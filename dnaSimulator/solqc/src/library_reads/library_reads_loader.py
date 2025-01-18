""""
This file is used to load the library file.
In case we want to enable different file types and format this loader should return a standardized Data Frame.
Should we want to process the file (drop columns, change upeer/lower etc..) these pre processing
events should happen here.
"""

import pandas as pd
from src.utils import biology
from src.utils.sequencing_file_parser import parse_sequencing_file


def get_filenames_from_file(file_name):
    with open(file_name, 'r') as file:
        files_names = file.read().splitlines()

    return files_names


def is_csv(value):
    if not isinstance(value, str):
        return False

    if value[-3:] != 'csv':
        return False

    return True


def retrieve_sequences(reads_input, n_reads, p):
    reads_input = get_filenames_from_file(reads_input)
    if isinstance(reads_input, list):
        sequences = []
        for file_name in reads_input:
            sequences.extend(parse_sequencing_file(file_name, n_reads, p))
            print(" - - Got sequences from {}".format(file_name))
    else:
        sequences = parse_sequencing_file(reads_input, n_reads, p)

    return sequences


def load_library(reads_input, as_df=True, n_reads=False, p=False):
    """

    :param reads_input: fastq files to get the sequences from or a csv of a pre matched library reads.
    :param as_df: If true return sequences as a pandas dataframe else return as a python list.
    :param n_reads: Optional value, will return the first n_reads if supplied.
    :param p: Optional value, if set will sample with rate p from the sequences.
    :return:
    """
    if is_csv(reads_input):
        return reads_input

    sequences = retrieve_sequences(reads_input, n_reads, p)

    if as_df:
        return pd.DataFrame(sequences, columns=['sequence'])
    else:
        return sequences
