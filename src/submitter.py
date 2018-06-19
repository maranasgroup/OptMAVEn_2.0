import cPickle as pkl
import os
import re
import subprocess
import sys
import tempfile

import benchmarking
import standards


class PbsBatchSubmitter(object):
    """ PbsBatchSubmitter creates and runs multiple jobs on the PBS queue.
    Jobs are grouped into 'batches' and each batch is run with one script.
    Thus, batches are run in parallel, while the jobs in one batch are run in series.
    """
    def __init__(self, experiment):
        # Set the attributes of the PbsBatchSubmitter.
        self.experiment = experiment
        self.purpose = experiment.purpose
        self.directory = tempfile.mkdtemp(dir=self.experiment.get_temp(), prefix="pbs_")
        self.file = os.path.join(self.directory, "pbs.pickle")
        self.jobs_file_prefix = os.path.join(self.directory, standards.PbsJobFilePrefix)
        self.jobs = dict()
        self.time_file_prefix = os.path.join(self.directory, "time-")
        self.array = 0
        self.callbacks = 0

    def save(self):
        """ Save the object as a pickle file. """
        with open(self.file, "w") as f:
            pkl.dump(self, f)

    def get_jobs_file(self, array):
        """ Get the name of the file in which are listed the jobs of the given job array. """
        return "{}{}".format(self.jobs_file_prefix, array)

    def get_time_file(self, array):
        """ Get the name of the file in which is written the time taken to run the given job array. """
        return "{}{}".format(self.time_file_prefix, array)

    def get_jobs_files(self):
        """ Get all of the job files. """
        if self.array > 0:
            return map(self.get_jobs_file, range(1, self.array + 1))
        else:
            return list()

    def get_time_files(self):
        """ Get all of the time files. """
        if self.array > 0:
            return map(self.get_time_file, range(1, self.array + 1))
        else:
            return list()

    def collect_times(self):
        """ Move the information from the time files into the Experiment's benchmarking record. """
        # Loop through all of the benchmarking time files.
        for i, _file in enumerate(self.get_time_files()):
            try:
                # Try to read the time information from the file.
                times = benchmarking.parse_time_file(_file)
            except IOError:
                # If that fails, do nothing, because the timing information isn't crucial to the experiment.
                pass
            else:
                # Create a Time task from the time information.
                task = benchmarking.Time(self.purpose, times, "Array {}".format(i))
                # Add that Time to the experiment's benchmarking record.
                self.experiment.add_benchmark(task)

    def collect_garbage(self):
        """ Remove all of the temporary files associated with the submission. """
        if os.path.isdir(self.directory):
            if self.array > 0 and self.experiment.benchmarking:
                # Collect benchmarking information before deleting the files.
                self.collect_times()
            self.experiment.safe_rmtree(self.directory)

    def submit(self, program, args, jobs):
        """ Submit the jobs as batches. """
        # jobs is a dictionary where the keys are the arguments that must be passed to experiment.py and the values are the files that are generated when the jobs have completed.
        self.program = program
        self.args = args
        self.jobs.update(jobs)
        # Remove all previous jobs files before submitting new jobs.
        self.collect_garbage()
        # We know a job is unfinished if the file that should be generated upon the job's completion does not exist.
        unfinished = [job for job, _file in self.jobs.iteritems() if not os.path.isfile(_file)]
        if len(unfinished) > 0:
            # If there are any unfinished jobs, make a directory in which the submitter can work.
            try:
                os.mkdir(self.directory)
            except OSError:
                pass
            # Count the number of batches (which are submitted as PBS job arrays).
            self.array = 0
            collection = list()
            # Loop through the unfinished jobs.
            for i, job in enumerate(unfinished):
                # Add the job to the list of jobs in the batch.
                collection.append(job)
                if len(collection) >= self.experiment.batch_size or i == len(unfinished) - 1:
                    # If the number of jobs in the batch equals the batch size, or if this is the last job:
                    # Increment the number of batches that must be run.
                    self.array += 1
                    # Name the file in which the names of the jobs in this batch are written.
                    jobs_file = self.get_jobs_file(self.array)
                    with open(jobs_file, "w") as f:
                        # Open the file and write a space-separated list of the names of the jobs in the batch.
                        f.write(" ".join(map(str, collection)))
                    # Erase the jobs from the collection so that the next batch starts empty.
                    collection = list()
            # Generate the command to run the jobs as a PBS job array.
            command = "{} {} {}{}".format(program, " ".join(map(str, args)), self.jobs_file_prefix, standards.PbsArrayId)
            if command.startswith("rm"):
                # To ensure that there are no rogue 'rm' processes, forbid submitting them.
                raise OSError("Submitting 'rm' is forbidden.")
            # Create a file for the script that will be submitted.
            handle, self.script_file = tempfile.mkstemp(dir=self.directory, prefix="script_", suffix=".sh")
            os.close(handle)
            if self.experiment.benchmarking:
                # If benchmarking is on, then re-format the command so that the run time is output.
                command = time_command(command, "{}{}".format(self.time_file_prefix, standards.PbsArrayId), standards.PbsTimeFormat)
            # Write the script.
            write_script(self.script_file, command, self.experiment.walltime, self.array)
            # Format the command to submit the script.
            call = [standards.PbsQsub, self.script_file]
            # Submit the script.
            p = subprocess.Popen(call, stdout=subprocess.PIPE)
            # Get the job array ID from the output of the submission.
            stdout, stderr = p.communicate()
            job_id_match = re.match("[0-9]+\[\]", stdout)
            if job_id_match is None:
                raise OSError("Unable to implement callback for job id {}".format(stdout))
            job_id = job_id_match.group()
            # Make a callback script that will run once the script that was just submitted has finished.
            # This script will run submitter.py and pass the submitter object's pickle file as an argument.
            # That command will cause this function (submit) to run again (see if __name__ == "__main__" at the bottom).
            # submit will check if all jobs have completed, and if so, initiate the next step in the experiment.
            self.callbacks += 1
            command = "{} {} {}".format(standards.PythonCommand, os.path.realpath(__file__), self.file)
            handle, self.callback_file = tempfile.mkstemp(dir=self.directory, prefix="callback_", suffix=".sh")
            os.close(handle)
            if self.experiment.benchmarking:
                # If benchmarking is on, then make a file to hold the information about how long the callback script ran.
                handle, self.callback_time_file = tempfile.mkstemp(dir=self.experiment.get_temp(), prefix="callback_time_", suffix=".txt")
                os.close(handle)
                self.experiment.add_time_file(self.callback_time_file, status_offset=1)
            else:
                # Otherwise, do not make such a file.
                self.callback_time_file = None
            # Save the submitter object.
            self.save()
            # Submit the callback script, ensuring that it will delay starting until the job array finishes.
            submit(self.callback_file, command, self.experiment.walltime, options={"-W": "depend=afteranyarray:{}".format(job_id)}, time_file=self.callback_time_file)
        else:
            # If all jobs have finished, then run the next part of the experiment.
            self.experiment.run_next()


def script_initial(queue, walltime, array=0):
    """ Create the text needed to set up a script for submission to the PBS queue. """
    # Convert the walltime (given in seconds) to hours:minutes:seconds.
    secs = int(walltime % 60)
    walltime_mins = int((walltime - secs) / 60)
    mins = walltime_mins % 60
    walltime_hours = (walltime_mins - mins) / 60
    hours = walltime_hours
    walltime_text = "{:0>2}:{:0>2}:{:0>2}".format(hours, mins, secs)
    # Redirect output to /dev/null.
    destination = "/dev/null"
    # Write the lines.
    lines = ["#!/bin/sh",
             "#PBS -A {}".format(queue),
             "#PBS -l walltime={}".format(walltime_text),
             "#PBS -j oe",
             "#PBS -o {}".format(destination)]
    if array > 0:
        # If submitting a job array, this extra line is needed.
        lines.append("#PBS -t {}-{}".format(1, array))
    # Stop the script upon encoutering errors.
    lines.append("set -euo pipefail")
    return lines


def time_command(command, time_file, time_format):
    """ Run a command using Unix time, which not only runs the command but also outputs how long it takes. """
    return "{} -o {} -f {} {}".format(standards.TimeCommand, time_file, time_format, command)


def write_script(file_name, command, walltime, array=0):
    """ Write a complete script that is submitted to the queue and runs the given command. """
    # Write the beginning part that sets up the PBS queueing.
    lines = script_initial(standards.PbsQueue, walltime, array)
    if isinstance(command, list):
        # If multiple commands are given (as a list), add them all.
        lines.extend(command)
    else:
        # Otherwise, just add one command.
        lines.append(command)
    # Join the commands, write them to a script, and make the script executible.
    script = "\n".join(map(str, lines))
    with open(file_name, "w") as f:
        f.write(script)
    subprocess.call(["chmod", "u+x", file_name])


def submit(file_name, command, walltime, options=None, queue=True, purpose=None, time_file=None):
    """ Submit a PBS script for a stand-alone job (not a job array). """
    if time_file is not None:
        # If the time that the job takes to run should be saved, do so with a time command.
        command = time_command(command, time_file, standards.PbsTimeFormat)
    # Write the script.
    write_script(file_name, command, walltime)
    if options is None:
        # If no command-line options for running the script are given, make an empty options dict.
        options = dict()
    if queue:
        # The script can be submitted on the PBS queue.
        call = [standards.PbsQsub]
    else:
        # Alternatively, it can be run from the current shell.
        call = list()
    # Add any command-line options to the call.
    [call.extend([k, v]) for k, v in options.items()]
    # Add the name of the script to the call.
    call.append(file_name)
    # Run or submit the script.
    status = subprocess.call(call)


if __name__ == "__main__":
    """ If submitter.py is run directly, unpickle and submit a PbsBatchSubmitter object.
    usage:
    path/to/python path/to/submitter.py path/to/submitter_object.pickle
    """
    # Load the PbsBatchSubmitter.
    sub_file = sys.argv[1]
    with open(sub_file) as f:
        sub = pkl.load(f)
    # Submit the PbsBatchSubmitter. This function will check if all jobs have completed.
    # If not, missing jobs are repeated. If so, the next step of the experiment is run.
    sub.submit(sub.program, sub.args, dict())
