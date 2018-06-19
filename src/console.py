""" This module controls output to the console. """

import os

import standards


def disp(text):
    """ A wrapper for print. """
    print(wrap(text))


def clear():
    """ Clear the console. """
    os.system("clear")


def wrap(text, width=None):
    """ Wrap text to a specified number of columns. """
    # I no longer think wrapping text is necessary, and actually might be detrimental.
    # Thus, this function just returns the input text.
    if width is None:
        width = standards.ConsoleWidth
    return text


def disp_list(list_, sep=", ", max_display=standards.MaxListDisplay, get_string=False):
    """ Display a list of items with a given separator. """
    if max_display is None:
        max_display = len(list_)
    max_display = min(max_display, len(list_))
    # Allow an upper limit on the number of items to print.
    if max_display < len(list_):
        items = list_[: max_display - 1] + ["...", list_[-1]]
    else:
        items = list_
    joined = sep.join(map(str, items))
    # The output can either be returned as a string (get_string is True) or printed directly (get_string is False).
    if get_string:
        return wrap(joined)
    else:
        disp(joined)
