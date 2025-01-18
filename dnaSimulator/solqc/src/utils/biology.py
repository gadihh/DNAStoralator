import re
import numpy as np
EQUAL = "0"
DELETION = "1"
INSERTION = "2"
MISMATCH = "3"

CIGAR_DICT = {
    "=" : "0",
    "I" : INSERTION,
    "D" : DELETION,
    "X" : MISMATCH
}


def get_reverse_compliment(sequence):
    bases_compliment = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    reverse = sequence[::-1]
    reverse_compliment = [bases_compliment[x] for x in reverse]
    return "".join(reverse_compliment)



def sample_fastq_sequences(path_to_file, p):
    sequences = []
    with open(path_to_file) as fastq_file:
        line = fastq_file.readline()
        while line:
            line = fastq_file.readline()
            if re.match("[ACGT]+\Z", line[:-1]) and np.random.binomial(1, p):
                sequences.append(line[:-1])

    return sequences

def n_fastq_sequences(path_to_file, n):
    sequences = []
    with open(path_to_file) as fastq_file:
        line = fastq_file.readline()
        while line:
            line = fastq_file.readline()
            if re.match("[ACGT]+\Z", line[:-1]):
                sequences.append(line[:-1])
                if n != -1 and len(sequences) >= n:
                    return sequences

    return sequences

def get_fastq_sequnces(path_to_file, n=False, p=False):
    if p:
        return sample_fastq_sequences(path_to_file, p)

    elif n:
        return n_fastq_sequences(path_to_file, n)

    else:
        return n_fastq_sequences(path_to_file, -1)


def parse_cigar(path):
    values = re.findall(r'\d+', path)
    chars = re.findall(r'[=IDX]', path)

    str = ""
    for val, char in zip(values, chars):
        tmp = CIGAR_DICT[char]*int(val)
        str += tmp

    return str


def str_cigar_path_recounstruct(s, path):
    reconstruct = s

    for i, char in enumerate(path):
        if char == EQUAL:
            pass

        if char == DELETION:
            reconstruct = reconstruct[:i] + reconstruct[i + 1:]

        if char == INSERTION:
            reconstruct = reconstruct[:i-1] + 'I' + reconstruct[i-1:]

        if char == MISMATCH:
            reconstruct = reconstruct[:i-1] + 'M' + reconstruct[i:]

    return reconstruct
