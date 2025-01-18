""" Configuration file.
This module serves as a global configuration file, which
is accessible horizontally across the program.
"""

import json

configuration_file = {}
deliverable_path = "deliverable/"

def set_configuration_file(file):
    """
    Set the user configuration file.
    :param file: The file which will be used as the configuration file.
    """
    global configuration_file
    configuration_file = file

def get_configuration():
    """
    Get the configuration file attributes as a python dict.
    :return: A dict containing all the attributes specified in the config file.
    """
    global configuration_file
    with open(configuration_file) as f:
        configuration = json.load(f)

    return configuration


def set_deliverable_path(time=""):
    global deliverable_path
    deliverable_path = "deliverable/{}".format(time).replace(":", "-")


def get_deliverable_path():
    return deliverable_path
