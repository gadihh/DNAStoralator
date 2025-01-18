'''
This model is in charge of parsing sequencing files.
'''
import numpy as np
import re

from enum import Enum
class SequencingFile(Enum):
    FASTQ = 1
    FASTA = 2


'''
**********************************************************************************************************************'
                                                    API
**********************************************************************************************************************'
'''


def parse_sequencing_file(sequencing_file_name, n=False, p=False):
    '''
    Given a sequencing file return a list of all the sequences in the file.
    Currently supports fastq and fastq.
    :param file: The name of the files holding the sequencing
    :return: list: a list containing all the sequences in the file.
    '''
    file_type = get_file_type(sequencing_file_name)
    if file_type == SequencingFile.FASTQ:
        return parse_fastq_file(sequencing_file_name, n, p)


def get_file_type(file_name):
    return SequencingFile.FASTQ


'''
**********************************************************************************************************************'
                                                    FastQ parsing
**********************************************************************************************************************'
'''
def get_fastq_sequences(file_name):
    sequences = []
    with open(file_name) as fastq_file:
        line = fastq_file.readline()
        while line:
            line = fastq_file.readline()
            if re.match("[ACGT]+\Z", line[:-1]):
                sequences.append(line[:-1])

    return sequences


def parse_fastq_file(file_name, n=False, p=False):
    sequences = get_fastq_sequences(file_name)
    if p:
        n = int (p * len(sequences))
        return np.random.choice(sequences, n, replace=False)

    elif n:
        return np.random.choice(sequences, n, replace=False)

    else:
        return sequences
'''
**********************************************************************************************************************'
                                                    FastA parsing
**********************************************************************************************************************'
'''

def read_fasta(fp):
    sequences = []
    seq = []
    name=None

    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                sequences.append(''.join(seq))
            seq = []
            name = line
        else:
            seq.append(line)

    if name:
        sequences.append(''.join(seq))

    return sequences


def get_fasta_sequences(file_name):
    with open(file_name) as fp:
        return read_fasta(fp)


def parse_fasta_file(file_name, n=False, p=False):
    sequences = get_fasta_sequences(file_name)

    if p:
        n = int(p * len(sequences))
        return np.random.choice(sequences, n, replace=False)

    elif n:
        return np.random.choice(sequences, n, replace=False)

    else:
        return sequences
