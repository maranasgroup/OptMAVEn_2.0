""" This module provides functions for user input. """

import copy
import os
import string

from Bio.PDB import PDBList
import numpy as np

from console import clear, disp, disp_list, wrap
import standards


def get(prompt):
    """ Get a string of input. """
    return raw_input(wrap(prompt))

def get_number(prompt, min_value=None, max_value=None, value_type=float):
    """ Get numeric input. Optionally have a min_value and max_value for the number, or force the number to be value_type = int. """
    if value_type not in (float, int):
        raise ValueError("Numeric type must be float or int, not {}.".format(value_type))
    do = True
    while do:
        try:
            answer = value_type(get(prompt))
        except ValueError:
            pass
        else:
            if min_value is not None and answer < min_value:
                disp("Value must be greater than or equal to {}.".format(min_value))
            elif max_value is not None and answer > max_value:
                disp("Value must be less than or equal to {}.".format(max_value))
            else:
                do = False
    return answer


def get_range(prompt, min_value=None, max_value=None, value_type=float):
    """ Get a range (np.array) of values. Optionally constrain the range between min_value and max_value, or force a value_type for the first and last numbers. """
    disp(prompt)
    first = get_number("First number: ", min_value, max_value, value_type)
    last = get_number("Last number: ", min_value, max_value, value_type)
    if first == last:
        return [first]
    values = get_number("Number of values: ", 1, None, int)
    return np.linspace(first, last, values).tolist()
    

def get_file(prompt, directory, new_file=False, fetch_pdb=False):
    """ Get a filename that exists within a given directory.
    new_file: if True, file must be a valid but non-existant path. 
    fetch_pdb: if True, if the file does not exist but is a valid PDB id, then download the PDB to directory. """
    # The directory must exist.
    if not os.path.isdir(directory):
        raise OSError("Directory does not exist: {}".format(directory))
    if new_file:
        # In this mode, a new file must be specified.
        do = True
        while do:
            base_file = get(prompt)  # name of the new file
            full_file = os.path.join(directory, base_file)  # full path to the new file
            if os.path.exists(full_file):
                # A new file may not overwrite an existing file.
                disp("Cannot create file, path already exists: {}".format(full_file))
            else:
                # Ensure the name of the file is a valid path.
                if standards.is_path_component(base_file):
                    file_directory = os.path.dirname(full_file)
                    # Ensure the directory in which the file will be written exists.
                    if os.path.isdir(file_directory):
                        # If so, then the file path is valid: return it.
                        do = False
                    else:
                        disp("Cannot create file, directory does not exist: {}".format(file_directory))
                else:
                    disp("Only the following characters are permitted in path components: {}".format(standards.AllowedPathCharacters))
        return full_file
    else:
        # In this mode, an existing file must be specified.
        do = True
        while do:
            name = get(prompt)  # name of the file
            if name in os.listdir(directory):
                # If the file exists, then nothing more must be done: return the file path.
                do = False
            elif fetch_pdb:
                try:
                    # If the file does not exist but fetch_pdb is True, then try to download the PDB.
                    PDBList().retrieve_pdb_file(name, pdir=directory, file_format="pdb")
                    # PDBs are downloaded in this file name format.
                    temp_name = os.path.join(directory, "pdb{}.ent".format(name))
                    # Rename the file to PDBID.pdb
                    name += ".pdb"
                    os.rename(temp_name, os.path.join(directory, name))
                    # If that worked, return the path to the PDB file.
                    do = False
                except:
                    disp("Could not fetch PDB '{}'".format(name))
            else:
                disp("File does not exist '{}'".format(name))
        # Return the path to the file.
        return os.path.join(directory, name)


def get_files(prompt, directory, min_number, max_number):
    """ Select between min_number and max_number of files from a list. """ 
    return [os.path.join(directory, file_) for file_ in select_from_list(prompt, os.listdir(directory), min_number, max_number)]


def get_yn(prompt):
    """ Get a yes/no response. Return True for yes, False for no. """
    yes = ["yes", "y"]
    no = ["no", "n"]
    answer = None
    while answer not in yes + no:
        answer = get(prompt).lower()
    return answer in yes


def parse_numeric_selection(selection_text):
    """ Return a set of whole numbers based on a selection string. The string should be formatted with commas and dashes, e.g. "1,3,7-9" returns [1, 3, 7, 8, 9]. Spaces do not matter. """
    # Split the selection text into items by ','
    selection_items = "".join(selection_text.split()).split(",")
    selection_numbers = list()
    # Loop through the items.
    for item in selection_items:
        if "-" in item:
            # A '-' in the item indicates a range.
            bounds = item.split("-")
            if len(bounds) != 2:
                # There must be exactly two numbers: one on each side of the '-'.
                raise ValueError('Cannot parse range: {}. Must be formatted "x-y".'.format(item))
            # Add all numbers in the range (inclusive) between the two given numbers.
            selection_numbers.extend(range(int(bounds[0]), int(bounds[1]) + 1))
        else:
            # Otherwise, just convert to an integer and append the item.
            selection_numbers.append(int(item))
    return selection_numbers


def select_from_list(prompt, list_, min_number=None, max_number=None, names=None):
    """ Select between min_number and max_number of items from a list. If names is given, list element names are replaced with names. """
    n = len(list_)
    if isinstance(list_, dict):
        # If the list is given as a dict:
        if names is None:
            # If names is not given, use the keys as names and the values as the items in the list.
	    names = map(str, list_.keys())
            list_ = list_.values()
        else:
            # Otherwise, discard the values and use the keys.
            list_ = list_.keys()
    if names is None:
        # If names is not specified, then convert the list items to strings and use them as names.
        names = map(str, list_)
    else:
        # Ensure that the number of names matches the number of items in the list.
        if len(names) != n:
            raise ValueError("Length of names must equal length of list.")
    # Strip off whitespace from the names.
    names = map(string.strip, names)
    # The minimum number selected cannot be negative.
    if min_number is None:
        min_number = 0
    else:
        min_number = max(int(min_number), 0)
    # The maximum number selected cannot exceed the total number of items.
    if max_number is None:
        max_number = n
    else:
        max_number = min(int(max_number), n)
    # The minimum number cannot exceed the maximum number.
    if min_number > max_number:
        raise IOError("Minimum selction {} cannot be greater than maximum selection {}.".format(min_number, max_number))
    # If the minumum number equals the number of items, all items must be selected.
    if min_number == n:
        disp("Automatically selecting all {n} items from length-{n} list: {l}".format(n=n, l=", ".join(names)))
        # Return a copy of the list so that things that depend on the lists input and returned being a different objects won't fail.
        return list(list_)
    # Initialize to no items selected.
    items = None
    # Make a statement of which options are available.
    permitted = "Options include the following:\n{}".format(disp_list(names, get_string=True))
    # Has that statment been displayed yet?
    permitted_displayed = False
    # Repeat until a selection has been made.
    while items is None:
        # Ask the user to specify items.
        selection_string = get(prompt)
        if selection_string.strip() == "":
            # If the input is blank, say what options are available.
            disp(permitted)
            permitted_displayed = True
        else:
            # Otherwise, try to parse the input.
            selection_keys = split_selection_string(selection_string)
            try:
                indexes = [index for key in selection_keys for index in get_selected_indexes(names, key)]
            except ValueError:
                # If parsing fails, then stop parsing and get another input from the user.
                select_items = False
            else:
                # If parsing succeeded, then proceed to the selection step.
                select_items = True
            if select_items:
                # Select the items.
                selected_items = [list_[index] for index in indexes]
                # Ensure the number selected meets specifications.
                if len(selected_items) >= min_number and len(selected_items) <= max_number:
                    # If so, then finalize the selection by setting items = selected_items, which terminates the loop.
                    items = selected_items
                else:
                    disp("You must select from {} to {} items.".format(min_number, max_number))
            else:
                # Alert that something went wrong with parsing.
                disp('There was a problem selecting the following items: "{}"'.format(selection_string))
                if not permitted_displayed:
                    # If an error occurs and the options have not been displayed yet, then display them.
                    disp(permitted)
                    permitted_displayed = True
    return items


def split_selection_string(string, sep=","):
    """ Split a string, then strip whitespace from each split item. """
    return [item.strip() for item in string.split(sep)]


def get_selected_indexes(list_, key, range_sep="-"):
    """ Return the indexes of the items specified in the selection syntax. """
    if key == standards.SelectAll:
        # If the keyword 'all' is used, then return all indexes.
        return range(len(list_))
    elif key == standards.SelectNone:
        # If the keyword 'none' is used, then return no indexes.
        return list()
    elif key.count(range_sep) == 1:
        # If the selection key indicates a range, then use the given items as bounds.
        lower, upper = [unescape_item(item.strip()) for item in key.split(range_sep)]
        # Return all items between them.
        i_lower = list_.index(lower)
        i_upper = list_.index(upper)
        return range(i_lower, i_upper + 1) 
    else:
        # Otherwise, just use the one index.
        return [list_.index(unescape_item(key))]


def unescape_item(item):
    """ Remove an escape character if an item starts with the escape character. """
    if item.startswith(standards.EscapeCharacter):
        return item[1:]
    else:
        return item


def select_one_from_list(prompt, list_, names=None):
    """ Select exactly one item from a list, and return that item, not a length-one list of that item. """
    return select_from_list(prompt, list_, 1, 1, names)[0]
