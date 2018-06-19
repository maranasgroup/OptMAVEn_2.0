""" Benchmarking provides functions for measuring the time and disk usage of an experiment. """


from collections import OrderedDict
import datetime
import os
import subprocess

import standards


class BlankTimeFileError(IOError):
    """ Raised whenever a file that should contain time information is blank. """
    pass


class Task(object):
    """ A Task is any process that an Experiment performs. The runtime or drive usage of a Task can be measured. """
    def __init__(self, purpose, detail=None, time_stamp=None):
        # purpose is a string that briefly describes the purpose of the Task.
        # For example, an interaction energy calculation step could have purpose = 'MAPs energy calculations'
        self.purpose = purpose
        # detail describes the Task more specifically than purpose.
        # For example, if purpose = 'MAPs energy calculations', a potential value for detail would be 'Array 1', indicating that the Task is about completing job array 1 of the MAPs energy calculations.
        self.detail = detail
        # The Task may be given a time stamp to record when it was created.
        # The time stamp does not affect the recorded duration of the Task, but it does tell when the Task began.
        if time_stamp is None:
            # If no time stamp is given, the time at which the Task was initialized is used.
            self.time_stamp = datetime.datetime.now()
        else:
            self.time_stamp = time_stamp

    def get_time_stamp(self):
        return self.time_stamp


class Time(Task):
    """ Time is a kind of Task whose duration is measured.
    Time objects are used for processes that are submitted as batch jobs (using submitter.PbsBatchSubmitter).
    Time objects are also created temporarily by TimeFile.to_dict. """
    def __init__(self, purpose, time, detail=None, time_stamp=None):
        Task.__init__(self, purpose, detail, time_stamp)
        # Time objects also include an attribute time, which records the amount of time the Task took.
        self.time = time

    def to_dict(self):
        # Convert the Time object to a dictionary containing the object's information in a readily usable format.
        info = {
            "Type": "Time",
            "Purpose": self.purpose,
            "Detail": self.detail,
            "Time Stamp": self.time_stamp.strftime(standards.DatetimeFormat)
        }
        # For each time code (real, sys, user), store the time taken.
        for _type in standards.UnixTimeCodes:
            info[_type] = self.time[_type]
        return info
            

class TimeFile(Task):
    """ TimeFile is a Task that measures the duration of a process.
    TimeFiles are used for jobs that are not submitted as batches (e.g. the initial script of an OptmavenExperiment).
    TimeFiles are associated with a file path, in which the time information is recorded.
    The time information is not present upon the creation of a TimeFile, but is rather added later when the timed process ends.
    """
    def __init__(self, _file, purpose, detail=None):
        Task.__init__(self, purpose, detail)
        # TimeFiles record the path to the file containing the time information.
        self.file = _file

    def to_dict(self):
        # Make a Time object from the contents of the file, then call to_dict() on the Time object.
        return Time(self.purpose, parse_time_file(self.file), self.detail, self.time_stamp).to_dict()
            

class DriveUsage(Task):
    """ DriveUsage Tasks record the drive usage (storage space) consumed by the experiment's directory at the time when the Task is created.
    """
    def __init__(self, experiment, detail=None):
        Task.__init__(self, experiment.purpose, detail)
        # Run the Unix command du to get the drive usage.
        du = subprocess.Popen(["du", "-s", experiment.directory], stdout=subprocess.PIPE)
        # There is no need to use an intermediate file because the output of du is obtained with subprocess.communicate.
        out, error = du.communicate()
        drive_usage, directory = out.split("\t")
        # Record the drive usage as a float.
        self.drive_usage = float(drive_usage)
        
    def to_dict(self):
        # Return a dict of the important information of the DriveUsage object.
        info = {
            "Type": "Drive Usage",
            "Purpose": self.purpose,
            "Detail": self.detail,
            "Drive Usage": self.drive_usage,
            "Time Stamp": self.time_stamp.strftime(standards.DatetimeFormat)
        }
        return info
 

def parse_time_file(_file):
    """ Open a file containing time information, validate the contents, and return a dict of the time codes and their values. """
    with open(_file) as f:
        lines = f.readlines()
    # Raise an error on blank files.
    if len(lines) == 0:
        raise BlankTimeFileError("Found blank time file: {}".format(_file))
    # Raise an error if the number of lines doesn't match the number of time codes, because there should be one time code per line.
    if len(lines) != len(standards.UnixTimeCodes):
        raise ValueError("Need {} lines in time file {}".format(len(standards.UnixTimeCodes), _file))
    # Return a dict of {time code: time value}.
    return {time_type: float(line) for line, time_type in zip(lines, standards.UnixTimeCodes)}
