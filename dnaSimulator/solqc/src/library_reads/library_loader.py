""""
This file is used to load the library file.
In case we want to enable different file types and format this loader should return a standardized Data Frame.
Should we want to process the file (drop columns, change upeer/lower etc..) these pre processing
events should happen here.
"""

import pandas as pd
from src.utils import biology
from src.utils.sequencing_file_parser import parse_sequencing_file


def load_library(file, as_df=True, n_reads=False, p=False):
    """

    :param file: A file or files of fastq to get the sequences from.
    :param as_df: If true return sequences as a pandas dataframe else return as a python list.
    :param n_reads: Optional value, will return the first n_reads if supplied.
    :param p: Optional value, if set will sample with rate p from the sequences.
    :return:
    """
    if isinstance(file, list):
        sequences = []
        for file_name in file:
            #sequences.extend(parse_sequencing_file(file_name, n_reads, p))
            sequences.extend(biology.get_fastq_sequnces(file_name, n_reads, p))
            print(" - - Got sequences from {}".format(file_name))
    else:
        #sequences = parse_sequencing_file(file, n_reads, p)
        sequences = biology.get_fastq_sequnces(file, n_reads, p)

    if as_df:
        return pd.DataFrame(sequences, columns=['sequence'])
    else:
        return sequences