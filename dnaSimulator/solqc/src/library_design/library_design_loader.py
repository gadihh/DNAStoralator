""""
This file is used to load the oligo design file.
In case we want to enable different file types and format this loader should return a standardized Data Frame.
Should we want to process the file (drop columns, change upeer/lower etc..) these pre processing
events should happen here.
"""

# Imports
import sys
import pandas as pd
from pathlib import Path

import src.config as config
import src.utils.error_logger as error_logger

def load_oligo_design(design_input):
    """
    Creates a dataframe containing the data relevant for the Oligo Design object.
    :param design_input: (Str) Either a file name pointing to a design csv or iupac string.
    :return A datafrmae containg either the data from the file or data derived from the iupac string.
    """
    exists = Path(design_input).is_file()
    if exists:
        return load_oligo_csv_design(design_input)
    else:
        return load_oligo_uipac_design(design_input)

def load_oligo_csv_design(file):
    """
    Reads the given file name into a pandas dataframe, returns the dataframe.
    :param file: The file name from which to build the dataframe.
    :return: A Dataframe containing the data in the file.
    """
    oligo_df = pd.read_csv(file)
    return oligo_df


def get_first_occurence_index(s, options):
    for i, letter in enumerate(s):
        if letter in options:
            return i


def get_last_occurence_index(s, options):
    for i, letter in enumerate(reversed(s)):
        if letter in options:
            return len(s) - i


def is_valid_word(s, vocab):
    for letter in s:
        if letter not in vocab:
            return False

    return True
def load_oligo_uipac_design(iupac_string):
    IUPAC = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': 'AG',
        'Y': 'CT',
        'S': 'GC',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }
    IUPAC_LETTERS = ['R','Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']

    if not is_valid_word(iupac_string, IUPAC.keys()):
        error_string = "There is a problem with the design input either supply a valid csv file or an iupac string."
        logger = error_logger.get_logger()
        logger.log_error(error_string)

        exit(1)

    # List to hold the different variants. Once the recuression (extracto code) will run this
    # List will hold all the possibel sequences.
    variants = []

    def char_array_to_string(c_array):
        return ''.join(c_array)

    def extract_code(word, i):
        if i == len(word):
            variants.append(char_array_to_string(word))
            return

        for letter in IUPAC[word[i]]:
            word[i] = letter
            extract_code(word.copy(), i + 1)

    extract_code(list(iupac_string), 0)
    df = pd.DataFrame(variants, columns=['sequence'])

    # iupac_start = get_first_occurence_index(iupac_string, IUPAC_LETTERS)
    # iupac_end = get_last_occurence_index(iupac_string, IUPAC_LETTERS)

    configuration = config.get_configuration()
    iupac_start = configuration['barcode_start']
    iupac_end = configuration['barcode_end']

    df['barcode'] = [sequence[iupac_start:iupac_end] for sequence in variants]

    return df

