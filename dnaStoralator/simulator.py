# Module: simulator.py
from PyQt5.QtCore import pyqtSignal

import strand_error_sim
import copy
import platform
import subprocess
import os
from scipy.stats import skewnorm
import edlib
from custom_random_variable import CustomRvContinuous
import random

class Simulator:
    """
    Simulator class.
    Holds the attributes needed for error simulation on the strands input file:

    Class variables:
    :var self.total_error_rates: Dictionary of the total error rates used in the simulation, as provided in error_rates
        parameter.
    :var self.base_error_rates: Error rates corresponding to each base, as passed.
    :var self.long_deletion_length_rates: Based on (excluding single base deletion):
        https://www.biorxiv.org/content/biorxiv/early/2019/11/13/840231/F12.large.jpg?width=800&height=600&carousel=1
    :var self.input_path: Path of input file, as passed.
    :var self.is_stutter_method: Indicates whether simulator should use stutter method or not.
    :var self.distribution_info: Dictionary containing the information of the desired distribution of the number of
        copies to generate for each strand. Initiated based on user input:
        - User can supply a dictionary (as explained below) with the the distribution type, function or list of values,
            minimum and maximum values.
            Possible dictionary values:
            + 'type' : 'continuous' for a random variable, 'vector' for a custom vector of values
            + 'value' : the string of the random variable (IN PYTHON SYNTAX) OR the vector of the custom values
            + 'min' : minimum number of copies to generate. required only for type 'continuous'.
            + 'max' : maximum number of copies to generate. required only for type 'continuous'.
        - Use default parameter to use default settings (normal distribution with skewness).
        Examples:
        * Random variable:
            {   'type': 'continuous',
                'value': 'x ** 2',  # Note the python syntax
                'min': 3,
                'max': 100
            }
        * Vector:
            {   'type': 'vector',
                'value': [1, 10, 5, 9, 6, 200]
            }
        * default:
            Pass nothing or None.

    :var self.random: Vector of number of copies for each origin strand. Initiated based on user input:
        - User can supply a vector of integers of the length of the input file (=number of origin strands in input file)
            which indicates how many copies will be generated for each strand. For example:
            If the input file has 5 origin strands, the vector can be: [1, 10, 3, 5, 2]
            (the 1st strand in the input will get 1 copy,
            the 2nd strand in the input will get 10 copies...)
        - User can supply a random variable. Use default parameter to use default settings (normal distribution with
            skewness).
    """
    def __init__(self, total_error_rates, base_error_rates, input_path, is_stutter_method=False, distribution_info=None):
        """
        @param total_error_rates: Dictionary of the total error rates used in the simulation.
            Example of a dictionary:
            {'d': 0.1, 'i': 0.2, 's': 0.1, 'ld': 0.6}
            NOTE: Dictionary can be passed with values as strings, as the constructor can to parse them to floats.
        @param base_error_rates: Dictionary of dictionaries for each base.
            Example:
            {   'A': {'s': 0.1, 'i': 0.2, 'pi': 0.1, 'd': 0.05, 'ld': 0.6},
                'T': {...},
                'C': {...},
                'G': {...}
            }
            NOTE: Dictionary can be passed with rates values as strings, as the constructor can to parse them to floats.
        @param input_path: path of the input file.
            Example of input file content:
            TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCC\n
            TGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGT\n
            AAATTTGCAACCAGAAATTTGCAACCAGAATTCACTAGAGGACGCACGCTCTATTTCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCC\n
            TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCC\n
            TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCC\n
        @param is_stutter_method: False by default, set to True if stutter method should be used instead of other
            methods.
            Note: other methods have the same error simulation algorithm but with different rates, therefore
            no other statements are needed except the rates.
        @param distribution_info: Dictionary containing the information of the desired distribution of the number of
                copies to generate for each strand. Initiated based on user input:
                - User can supply a dictionary (as explained below) with the the distribution type, function or list of values,
                    minimum and maximum values.
                    Possible dictionary values:
                    + 'type' : 'continuous' for a random variable, 'vector' for a custom vector of values
                    + 'value' : the string of the random variable (IN PYTHON SYNTAX) OR the vector of the custom values
                    + 'min' : minimum number of copies to generate. required only for type 'continuous'.
                    + 'max' : maximum number of copies to generate. required only for type 'continuous'.
                - Use default parameter to use default settings (normal distribution with skewness).
                Examples:
                * Random variable:
                    {   'type': 'continuous',
                        'value': 'x ** 2',  # Note the python syntax
                        'min': 3,
                        'max': 100
                    }
                * Vector:
                    {   'type': 'vector',
                        'value': [1, 10, 5, 9, 6, 200]
                    }
                * default:
                    Pass nothing or None.
        """
        self.total_error_rates = copy.deepcopy(total_error_rates)
        self.base_error_rates = copy.deepcopy(base_error_rates)
        parse_rates_dictionary(self.total_error_rates)
        parse_rates_dictionary(self.base_error_rates)
        self.long_deletion_length_rates = {2: 2.8 * (10 ** (-4)),
                                           3: 7.75 * (10 ** (-5)),
                                           4: 3.25 * (10 ** (-5)),
                                           5: 10 ** (-6),
                                           6: 5.5 * (10 ** (-8))}
        self.input_path = input_path
        self.is_stutter_method = is_stutter_method

        self.distribution_info = distribution_info
        self.random = None
        self.min_copies = 1
        self.max_copies = 499

    def simulate_errors(self, report_func, evyat_path, shuffled_path):
        """
        Simulates strands duplication with errors on the strands from the input file.
        Writes the output in HeadEvyaLuis.txt file in the following format:
            [original strand][\n]
            *****************************[\n]
            [copy][\n]
            [copy][\n]
            ...
            [copy][\n]
            [\n]
            [\n]
            [original strand][\n]
            *****************************[\n]
            [copy][\n]
            [copy][\n]
            ...
            [copy][\n]

            ...
        """

        # count num of origin strands, each line is a separate strand:
        num_values = 0
        with open(self.input_path, 'r') as input_f:
            for line in input_f:
                num_values += 1

        # TODO: allow user input
        # TODO: allow user-defined function OR free random vector
        if self.distribution_info is None:
            # generate number of copies for each strand, as the number of strands:
            # https://stackoverflow.com/questions/24854965/create-random-numbers-with-left-skewed-probability-distribution
            self.max_copies = 499
            skewness = 10  # Negative values are left skewed, positive values are right skewed.
            self.random = skewnorm.rvs(a=skewness, loc=self.max_copies, size=num_values)  # Skewnorm function
            self.random = self.random - min(self.random)  # Shift the set so the minimum value is equal to zero.
            self.random = self.random / max(self.random)  # Standardize all the values between 0 and 1.
            self.random = self.random * self.max_copies  # Multiply the standardized values by the maximum value.
            self.random = self.random + self.min_copies  # avoid below minimum
            self.random = [round(x) for x in self.random]  # convert to integers

        elif self.distribution_info['type'] == 'continuous':
            raw_samples = CustomRvContinuous.multithreaded_rvs(size=num_values,
                                                               pdf_str=self.distribution_info['value'],
                                                               min_value=self.distribution_info['min'],
                                                               max_value=self.distribution_info['max'])
            self.random = [round(x) for x in raw_samples]

        elif self.distribution_info['type'] == 'vector':
            self.random = self.distribution_info['value']
            self.min_copies = min(self.random)
            self.max_copies = max(self.random)

        # for each strand, copy it the corresponding generated number of times and simulate error on each copy:
        i = 0
        with open(self.input_path, 'r') as input_f:
            os.makedirs('./output', exist_ok=True)  # required to create output directory if it doesn't exist
            with open(evyat_path, 'w', newline="\n") as output_f:
                for line in input_f:
                    report_func(num_values, i)
                    # write ORIGINAL strand with divider first:
                    original_strand = line
                    # strip the strand from newline:
                    original_strand = original_strand.rstrip()
                    output_f.write(original_strand + '\n' + '*****************************\n')

                    # set the number of copies for each design:
                    num_copies = self.min_copies
                    # in case of user defined vector, vector can be shorter of longer than the real number of designs.
                    # so in case of a shorter vector, use all the available values, and generate random ones between the minimum and maximum for the rest:
                    if i < len(self.random):
                        num_copies = self.random[i]
                    else:
                        num_copies = random.randint(self.min_copies, self.max_copies + 1)

                    # for each strand, do the simulation on a copy of it num_copies (the generated number of copies) times:
                    for j in range(num_copies):

                        # duplicate strand to create a working (output) strand:
                        output_strand = copy.deepcopy(original_strand)
                        # create a strand simulator for it:
                        strand_error_simulator = strand_error_sim.StrandErrorSimulation(self.total_error_rates,
                                                                                        self.base_error_rates,
                                                                                        self.long_deletion_length_rates,
                                                                                        output_strand)
                        # simulate according to method:
                        if self.is_stutter_method:
                            output_strand = strand_error_simulator.simulate_stutter_errors_on_strand()
                        else:
                            output_strand = strand_error_simulator.simulate_errors_on_strand()

                        output_f.write(output_strand + '\n')

                    # after each strand, add 2 newlines:
                    output_f.write('\n\n')

                    i += 1

        # mess the order of the output strands into a new file:
        mess_output_strands(evyat_path, shuffled_path)


def pseudo_cluster(start, end, dist):
    """
    Performs pseudo clustering on a given evyat.txt file. Modifies the file.
    :param start: start index of the cluster's identifier sequence.
    :param end: end index (not inclusive) of the cluster's identifier sequence.
    :param dist: maximal edit distance allowed for a sequence to be included in a cluster.
    """
    with open('output/evyat.txt', 'r') as evyat_f:
        with open('output/evyat_temp.txt', 'w', newline='\n') as temp_f:

            line = evyat_f.readline()
            design = line

            while line:
                # save designs and file structure and filter only copies according to the parameters:
                temp_f.write(line)

                if line == '*****************************\n':
                    design_barcode = design[start:end]
                    # read next line to get potential copy:
                    line = evyat_f.readline()
                    while line != '\n':
                        # if the line is a copy, strip it from newline:
                        copy_strand = line
                        # strip the strand from newline:
                        copy_strand = copy_strand.rstrip()
                        copy_strand_barcode = copy_strand[start:end]
                        edit_distance = edlib.align(design_barcode, copy_strand_barcode)['editDistance']

                        # save stand if edit distance is less than the maximal allowed (or drop it if not)
                        if edit_distance <= dist:
                            temp_f.write(line)

                        # read next potential copy:
                        line = evyat_f.readline()

                    temp_f.write(line)  # write the newline (it will be skipped just before next outer loop iteration)

                else:
                    design = line  # if next line is a separator, last line stored is the design

                # read next:
                line = evyat_f.readline()

    os.remove('output/evyat.txt')
    os.rename(r'output/evyat_temp.txt', r'output/evyat.txt')


def mess_output_strands(evyat_path, shuffled_path):
    """
    Messes the output strands.
    Creates a temporary file from the evyat.txt to run the shuffle program on it,
    and creates a new output file errors_shuffled.txt with all output strands shuffled and not clustered.
    """
    output_f = open('output/errors_temp.txt', 'w', newline='\n')
    with open(evyat_path, 'r') as errors_f:
        try:
            # skip first line and ****:
            next(errors_f)
            next(errors_f)

            # iterate over output strands and skip origin and **** each time.
            # at this point, current line is the first output strand.
            for line in errors_f:
                if line == '\n':
                    # skip current line and another \n + 2 lines:
                    next(errors_f)  # first \n skipped, pointing to another \n now
                    next(errors_f)  # second \n skipped, pointing to origin strand now
                    next(errors_f)  # origin strand skipped, pointing to **** line.
                    # now next iteration will point to next output strand
                else:
                    output_f.write(line)
        except StopIteration:
            output_f.close()
    output_f.close()

    if platform.system() == "Linux":
        # linux
        args = ['shuf', 'output/errors_temp.txt', '-o', shuffled_path]
        subprocess.run(args)
    elif platform.system() == "Darwin":
        # OS X
        args = ['./shuffle_prog/shuf_mac', 'output/errors_temp.txt', '-o', shuffled_path]
        subprocess.run(args)
    elif platform.system() == "Windows":
        args = ['./shuffle_prog/shuf_windows.exe', 'output/errors_temp.txt', '-o', shuffled_path]
        subprocess.run(args)
    os.remove('output/errors_temp.txt')


def parse_rate(rate_str) -> float:
    """
    Parses string of a single string of a numeric value representing a rate.
    If the parameter is already a float, does nothing (and returns it).
    :param rate_str: Can be either:
        - a number, eg. '0.053'
        or
        - use 10^x exp, eg. '53E-4' (which is equivalent to 0.053 or 53 * 10^(-4))
        NOTE: it can be a float as well, in this case it will simply be returned.
    :return: The float value of the given string.
        NOTE: If a float is passed as argument, it'll be returned as-is.
    """
    if isinstance(rate_str, float):
        return rate_str
    # if string is a number, convert it to float as-is.
    # if string is represented with E, convert it to a number first.
    index = rate_str.find('E')
    if index == -1:
        return float(rate_str)
    else:
        num = float(rate_str[:index])
        exp = 10 ** float(rate_str[index + 1:])
        return num * exp


def parse_rates_dictionary(rates_dict):
    """
    Parses dictionary's rates to floats if they appear as strings.
    Handles both "one level" and "two level" dictionaries used in the simulator.
    Modifies the passed dictionary!
    :param rates_dict: The dictionary to parse.
    """
    for key, value in rates_dict.items():
        if isinstance(value, dict):  # the given dictionary is a base error rates dictionary
            # key is the base, value is the dictionary for the base, consisting of errors & rates.
            for error, rate in value.items():
                rates_dict[key][error] = parse_rate(rate)
        else:  # the given dictionary is a one-level dictionary, assuming there are no more types of dictionaries.
            rates_dict[key] = parse_rate(value)


# Testing:
#
# if __name__ == '__main__':
#     pseudo_cluster(1, 10, 5)
#
#     error_rates_example = {'d': 9.58 * (10 ** (-4)),
#                            'ld': 2.33 * (10 ** (-4)),
#                            'i': 5.81 * (10 ** (-4)),
#                            's': 1.32 * (10 ** (-3))}
#     base_error_rates_example = {'A':
#                                 {'s': 0.135 * (10**(-2)),
#                                  'i': 0.057 * (10**(-2)),
#                                  'pi': 0.059 * (10**(-2)),
#                                  'd': 0.099 * (10**(-2)),
#                                  'ld': 0.024 * (10**(-2))},
#                                 'C':
#                                     {'s': 0.135 * (10 ** (-2)),
#                                      'i': 0.059 * (10 ** (-2)),
#                                      'pi': 0.058 * (10 ** (-2)),
#                                      'd': 0.098 * (10 ** (-2)),
#                                      'ld': 0.023 * (10 ** (-2))},
#                                 'T':
#                                     {'s': 0.126 * (10 ** (-2)),
#                                      'i': 0.059 * (10 ** (-2)),
#                                      'pi': 0.057 * (10 ** (-2)),
#                                      'd': 0.094 * (10 ** (-2)),
#                                      'ld': 0.023 * (10 ** (-2))},
#                                 'G':
#                                     {'s': 0.132 * (10 ** (-2)),
#                                      'i': 0.058 * (10 ** (-2)),
#                                      'pi': 0.058 * (10 ** (-2)),
#                                      'd': 0.096 * (10 ** (-2)),
#                                      'ld': 0.023 * (10 ** (-2))}}
#
#     input_path_example = 'input/strands_in.txt'
#
#     sim = Simulator(error_rates_example, base_error_rates_example, input_path_example)
#
#     sim.simulate_errors()
#
#     sim = Simulator(error_rates_example, base_error_rates_example, input_path_example, True)
#
#     sim.simulate_errors()
#
#     distribution_example = {'type': 'vector', 'value': [100, 10, 5, 9, 6, 200]}
#
#     sim = Simulator(error_rates_example, base_error_rates_example, input_path_example, False, distribution_example)
#
#     sim.simulate_errors()
#


