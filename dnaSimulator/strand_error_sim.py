# Module: strand_error_sim.py

import random
import copy

"""
How to use:

The simulator gets a strand and implements an error.

1. Create an object of the simulator with a dictionaries of errors and the source strand passed as arguments:
    example of an error dictionary: {'d': 0.1, 'i': 0.2, 's': 0.1, 'ld': 0.6}
    example of a base error dictionary:
        {   'A': {'s': 0.1, 'i': 0.2, 'pi': 0.1, 'd': 0.05, 'ld': 0.6},
            'T': {...},
            'C': {...},
            'G': {...}
        }
    example of s strand: "TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCC"
        (can be any length)
        NOTE: SUPPLY STRAND ONLY, WITHOUT ANY OTHER CHARACTERS (INCLUDING NEWLINES)!
                               _  _
                             _/0\/ \_
                    .-.   .-` \_/\0/ '-.
                   /:::\ / ,_________,  \
                  /\:::/ \  '. (:::/  `'-;
                  \ `-'`\ '._ `"'"'\__    \
                   `'-.  \   `)-=-=(  `,   |
                jgs    \  `-"`      `"-`   /
        https://www.asciiart.eu/television/sesame-street
        
2. For each strand of your input file, call simulate_error_on_strand.
    For stutter method, call simulate_stutter_errors_on_strand.
3. WALLAH!


Strings of types of errors supported (in dictionary):

'd'     - deletion (one base)
'ld'    - long deletion (multiple base deletion)
'i'     - insertion
's'     - substitution

For base error rates: Same as above, adding:
'pi'     - symbol pre-insertion

For stutter error rates:
use the 'd' (deletion), 'i' (insertion) and 's' (substitution) keys.

Do not use:
'n'     - none (no error), used in some cases to represent a "dummy" error. (to receive no-error)
"""


class StrandErrorSimulation:
    """
    # Strand Error Simulation class.
    Holds the attributes needed for error simulation on a single strand:

    Class variables:
    :var self.total_error_rates: dictionary of the total error rates used in the simulation, as provided in error_rates
        parameter.
    :var self.base_error_rates: error rates corresponding to each base, as passed.
    :var self.deletion_length_rates: as passed.
        https://www.biorxiv.org/content/biorxiv/early/2019/11/13/840231/F15.large.jpg?width=800&height=600&carousel=1
        https://www.biorxiv.org/content/biorxiv/early/2019/11/13/840231/F12.large.jpg?width=800&height=600&carousel=1
    :var self.strand: the strand to simulate the error on, as passed. This is also the final strand.
    :var self.index: the index to implement the error on.
          Initialized to 0.
    """
    def __init__(self, total_error_rates, base_error_rates, deletion_length_rates, strand):
        """
        :param total_error_rates: Dictionary of the total error rates used in the simulation.
            Example of a dictionary:
            {'d': 0.1, 'i': 0.2, 's': 0.1, 'ld': 0.6}
        :param base_error_rates: Dictionary of dictionaries for each base.
            Example:
            {   'A': {'s': 0.1, 'i': 0.2, 'pi': 0.1, 'd': 0.05, 'ld': 0.6},
                'T': {...},
                'C': {...},
                'G': {...}
            }
        :param deletion_length_rates: Dictionary of deletion length rates for higher lengths than 1 (start from 2):
            Example:
            {2: 0.2, 3: 0.1, 4: 0.3, 5: 0.05, 6: 0.001}
        :param strand: The strand to implement errors on.
        """
        self.total_error_rates = total_error_rates
        self.base_error_rates = base_error_rates
        self.deletion_length_rates = deletion_length_rates
        self.strand = strand
        self.index = 0
        # for testing only:
        self.err_type = None

    ''' Main Methods: '''

    def simulate_errors_on_strand(self) -> str:
        """
        Simulates errors on the given strand and returns the target strand.
        Use for any method EXCEPT stutter.
        :return:
        Modified strand after errors simulation
        """
        while self.index < len(self.strand):
            self.simulate_error_on_base()
            self.index += 1
        return self.strand

    def simulate_stutter_errors_on_strand(self) -> str:
        """
        Simulates stutter errors on the given strand and returns the target strand.
        Use only for stutter.
        :return:
        Modified strand after errors simulation
        """
        while self.index < len(self.strand):
            self.simulate_stutter_error_on_base()
            self.index += 1
        # after implementing the stutter errors, go over the strand again and inject substitutions according to rates:
        self.index = 0
        while self.index < len(self.strand):
            base = self.strand[self.index]
            base_substitution_rate = self.base_error_rates[base]['s']
            options = ['y', 'n']
            rates = [base_substitution_rate, 1 - base_substitution_rate]
            draw = random.choices(options, weights=rates, k=1)
            if draw[0] == 'y':
                self.strand = self.inject_substitution()
            self.index += 1
        return self.strand

    ''' Helper Methods: '''

    def simulate_stutter_error_on_base(self):
        """
        Simulates a stutter error on the current base (in the current index)
        Modifies the working strand stored in class.
        NO RETURN VALUE
        """
        # this method doesn't have long deletion
        base = self.strand[self.index]
        base_deletion_rate = self.base_error_rates[base]['d']
        base_stutter_rate = self.base_error_rates[base]['i']
        # draw whether there was deletion or not:
        options = ['y', 'n']
        rates = [base_deletion_rate, 1 - base_deletion_rate]
        draw = random.choices(options, weights=rates, k=1)
        if draw[0] == 'y':
            self.strand = self.inject_error('d')
            self.err_type = 'd'  # for testing
            return
        else:
            rates = [base_stutter_rate, 1 - base_stutter_rate]
            draw = random.choices(options, weights=rates, k=1)
            is_stutter = (draw[0] == 'y')
            # for testing, count stutter times
            self.err_type = 0
            while is_stutter:
                self.strand = self.strand[:self.index] + base + self.strand[self.index:]
                # for testing, count stutter times
                self.err_type += 1
                # draw again before next iteration:
                draw = random.choices(options, weights=rates, k=1)
                is_stutter = (draw[0] == 'y')
            # increment index to approach next original base:
            self.index += 1
            # for testing:
            if self.err_type == 0:
                self.err_type = 'n'
            return

    def simulate_error_on_base(self):
        """
        Simulates any error EXCEPT stutter on the current base (in the current index).
        Modifies the working strand stored in class.
        """
        base = self.strand[self.index]
        # 1. summarize all error rates into total rate, and conclude the complementary non-error rate:
        total_error_rate = 0
        for value in self.total_error_rates.values():
            total_error_rate += value
        no_error_rate = 1 - total_error_rate
        assert(no_error_rate >= 0)

        # 2. draw whether there's error or not in the given rates:
        # Give higher rates to first third
        options = ['y', 'n']
        rates = [total_error_rate, no_error_rate]
        if self.index <= (1/3) * len(self.strand):
            rates[0] = rates[0] * 3/2
            rates[1] = 1 - rates[0]
        else:
            rates[0] = rates[0] * 3/4
            rates[1] = 1 - rates[0]
        draw = random.choices(options, weights=rates, k=1)

        # 3. check type of drawn result:
        # 3.1. If there's error:
        if draw[0] == 'y':
            # generate an error type:
            error_type = self.generate_error_type_for_base(base)
            self.err_type = error_type  # for testing only
            self.strand = self.inject_error(error_type)
        # 3.2. If there's no error - do nothing.
        else:
            self.err_type = 'n'  # for testing only

    def generate_error_type_for_base(self, base) -> str:
        """
        Generate an error from the error rates dictionary passed as arguments:
        Returns the error type generated for the base (string).
        :param base: The value of the base currently working on: 'A', 'T', 'C', 'G'.
        :return 'd': for deletion
        :return 'ld': for long deletion
        :return 'pi': for insertion (base on pre-insertion symbol)
        :return 's': for substitution
        """
        # create two lists of the dictionary - options list and rates list:
        options = []
        rates = []
        base_rates = copy.deepcopy(self.base_error_rates[base])

        # remove insertion rates (as they are not needed in this stage - they are used for inserted base generation)
        del base_rates['i']

        for key, value in base_rates.items():
            options.append(key)
            rates.append(value)

        #  Note: choices uses weights, and thus is equivalent to conditional probability!
        draw = random.choices(options, weights=rates, k=1)
        # return error type string:
        return draw[0]

    def inject_deletion(self, error_type) -> str:
        """
        Inject deletion to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :param error_type: Type of deletion: 'd' for single base deletion or 'ld' for long deletion.
        :return: a strand with the injected deletion.
        """
        modified_strand = ""

        if error_type == 'd':
            # single base deletion:
            if self.index == len(self.strand) - 1:
                modified_strand = self.strand[:self.index]
            else:
                modified_strand = self.strand[:self.index] + self.strand[self.index + 1:]

        elif error_type == 'ld':
            # multiple base deletion:
            long_del_dict = copy.deepcopy(self.deletion_length_rates)

            # draw length:
            options = list(long_del_dict.keys())
            rates = list(long_del_dict.values())
            draw = random.choices(options, weights=rates, k=1)

            deletion_length = draw[0]
            if self.index + deletion_length > len(self.strand) - 1:
                modified_strand = self.strand[:self.index]
            else:
                modified_strand = self.strand[:self.index] + self.strand[self.index + deletion_length:]

        # keep index the same! The original base in that index (or further) was deleted.
        self.index -= 1
        return modified_strand

    def inject_insertion(self) -> str:
        """
        Inject insertion to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :return: a strand with the injected insertion.
        """
        base_insertion_rates = {'A': self.base_error_rates['A']['i'],
                                'T': self.base_error_rates['T']['i'],
                                'C': self.base_error_rates['C']['i'],
                                'G': self.base_error_rates['G']['i']}
        options = list(base_insertion_rates.keys())
        rates = list(base_insertion_rates.values())
        draw = random.choices(options, weights=rates, k=1)
        modified_strand = self.strand[:self.index] + draw[0] + self.strand[self.index:]
        # increment index to approach next original base:
        self.index += 1
        return modified_strand

    def inject_substitution(self) -> str:
        """
        Inject substitution to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :return: a strand with the injected substitution.
        """
        base = self.strand[self.index]
        modified_strand = list(self.strand)
        bases = ['A', 'T', 'G', 'C']
        options = []
        for b in bases:
            if b != base:
                options.append(b)
        # Note: 'options' is defined by 'bases' so the order is always the same as in 'bases'.
        # Set rates according to the base:
        rates = [1, 1, 1]
        # if modified_strand[index] == 'G':
        #     rates = []
        draw = random.choices(options, weights=rates, k=1)
        modified_strand[self.index] = draw[0]
        modified_strand = ''.join(modified_strand)
        return modified_strand

    def inject_error(self, error_type: str) -> str:
        """
        Inject the error type to the given strand, starting from the base in `index` location of the strand.
        Returns a strand with the injected error.
        :param error_type: Error type to inject ('d', 'ld', 's', 'pi' as documented)
        :return: a strand with the injected error.
        """
        # check error type and act accordingly:
        if error_type == 'd' or error_type == 'ld':
            return self.inject_deletion(error_type)
        elif error_type == 'pi':  # pre insertion rates are rates for insertion error
            return self.inject_insertion()
        elif error_type == 's':
            return self.inject_substitution()


# Testing:

# if __name__ == '__main__':
#
#     # stutter test:
#     stutter_error_rates_example = {}  # can be empty - not used!
#     stutter_deletion_length_rates_example = {}  # can be empty - not used!
#
#     # base error rates can have only relevant values:
#     # NOTE: THESE VALUES ARE NOT REAL (they don't represent correct rates regarding conditional probability and so),
#     # TO CHECK CORRECTNESS IT IS NEEDED TO HAVE REAL VALUES AND THEIR TOTAL PROBABILITIES SUM
#     stutter_base_error_rates_example = {'A':
#                                         {'d': 0.1,
#                                          'i': 0.3},
#                                         'C':
#                                             {'d': 0.1,
#                                              'i': 0.2},
#                                         'T':
#                                             {'d': 0.1,
#                                              'i': 0.05},
#                                         'G':
#                                             {'d': 0.1,
#                                              'i': 0.05}}
#
#     # http://www.faculty.ucr.edu/~mmaduro/random.htm
#     example_strand = "TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTA" \
#                      "GTGCAATGGGGCTTTTTTTCCATAGTCCTCGAGAGGAGGAGACGTCAGTCCAGATATCTTTGATGTCGTGATTGGAAGGACCCTTGGCCCTCCA" \
#                      "CCCTTAGGCAGT"
#
#     # full stutter test:
#
#     full_stutter_sim_f = open('full_stutter_simulation', 'w')
#     full_stutter_err_type_f = open('full_stutter_error_types', 'w')
#     full_stutter_loc_f = open('full_stutter_locations', 'w')
#     for j in range(1000):
#         stutter_simulator = StrandErrorSimulation(stutter_base_error_rates_example, stutter_base_error_rates_example,
#                                                   stutter_base_error_rates_example, example_strand)
#         while stutter_simulator.index < len(stutter_simulator.strand):
#             stutter_simulator.simulate_stutter_error_on_base()
#             type_result_to_write = str(stutter_simulator.err_type) + '\n'
#             full_stutter_err_type_f.write(type_result_to_write)
#             if stutter_simulator.err_type != 'n':
#                 location_result_to_write = str(stutter_simulator.index) + '\n'
#                 full_stutter_loc_f.write(location_result_to_write)
#             output_strand_to_write = str(stutter_simulator.strand) + '\n'
#             full_stutter_sim_f.write(output_strand_to_write)
#             stutter_simulator.index += 1
#         full_stutter_sim_f.write(stutter_simulator.strand + '\n')
#
#     full_stutter_sim_f.close()
#     full_stutter_loc_f.close()
#     full_stutter_err_type_f.close()
#
#     # analyze types:
#
#     full_stutter_err_type_f = open('full_stutter_error_types', 'r')
#     hist = [['d', 0], ['i', 0], ['n', 0]]
#     lines = full_stutter_err_type_f.readlines()
#     for line in lines:
#         if line == 'd\n':
#             hist[0][1] += 1
#         if line == 'n\n':
#             hist[2][1] += 1
#         else:  # number appears
#             hist[1][1] += 1
#     full_stutter_err_type_f.close()
#
#     full_stutter_err_type_ana_f = open('stutter_error_types_analysis', 'w')
#     full_stutter_err_type_ana_f.write('d appearance rate: ' + str(hist[0][1] / (1000 * len(example_strand))) + '\n')
#     full_stutter_err_type_ana_f.write('i appearance rate: ' + str(hist[2][1] / (1000 * len(example_strand))) + '\n')
#     full_stutter_err_type_ana_f.write('n appearance rate: ' + str(hist[1][1] / (1000 * len(example_strand))) + '\n')
#     full_stutter_err_type_ana_f.close()
#
#     # 'd': 0.0009580000000000001
#     # 'ld': 0.00023300000000000003
#     # 's': 0.00132
#     # 'i': 0.000581
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
#     deletion_length_rates_example = {2: 2.8 * (10 ** (-4)),
#                                      3: 7.75 * (10 ** (-5)),
#                                      4: 3.25 * (10 ** (-5)),
#                                      5: 10 ** (-6),
#                                      6: 5.5 * (10 ** (-8))}
#
#     simulator = StrandErrorSimulation(error_rates_example, base_error_rates_example, deletion_length_rates_example,
#                                       example_strand)
#
#     # error type & location simulation test:
#
#     err_type_f = open('error_types', 'w')
#     error_loc_f = open('error_locations', 'w')
#     output_strand_f = open('output_strand', 'w')
#     while simulator.index < len(simulator.strand):
#         simulator.simulate_error_on_base()
#         type_result_to_write = str(simulator.err_type) + '\n'
#         err_type_f.write(type_result_to_write)
#         if simulator.err_type != 'n':
#             location_result_to_write = str(simulator.index) + '\n'
#             error_loc_f.write(location_result_to_write)
#         output_strand_to_write = str(simulator.strand) + '\n'
#         output_strand_f.write(output_strand_to_write)
#         simulator.index += 1
#
#     err_type_f.close()
#     error_loc_f.close()
#     output_strand_f.close()
#
#     # General single strand test:
#
#     simulator = StrandErrorSimulation(error_rates_example, base_error_rates_example, deletion_length_rates_example,
#                                       example_strand)
#     full_strand_sim_f = open('full_strand_simulation', 'w')
#     full_strand_sim_f.write('original:\n' + example_strand + '\n')
#     result = simulator.simulate_errors_on_strand()
#     result_to_write = 'modified:\n' + result + '\n'
#     full_strand_sim_f.write(result_to_write)
#     full_strand_sim_f.close()
#
#     # full copy error types and locations test:
#
#     full_strand_copy_sim_f = open('full_strand_copy_simulation', 'w')
#     full_copy_err_type_f = open('full_copy_error_types', 'w')
#     full_copy_error_loc_f = open('full_copy_error_locations', 'w')
#     for j in range(1000):
#         simulator = StrandErrorSimulation(error_rates_example, base_error_rates_example, deletion_length_rates_example,
#                                           example_strand)
#         while simulator.index < len(simulator.strand):
#             simulator.simulate_error_on_base()
#             type_result_to_write = str(simulator.err_type) + '\n'
#             full_copy_err_type_f.write(type_result_to_write)
#             if simulator.err_type != 'n':
#                 location_result_to_write = str(simulator.index) + '\n'
#                 full_copy_error_loc_f.write(location_result_to_write)
#             output_strand_to_write = str(simulator.strand) + '\n'
#             full_strand_copy_sim_f.write(output_strand_to_write)
#             simulator.index += 1
#         full_strand_copy_sim_f.write(simulator.strand + '\n')
#
#     full_copy_err_type_f.close()
#     full_copy_error_loc_f.close()
#     full_strand_copy_sim_f.close()
#
#     # analyze types:
#
#     full_copy_err_type_f = open('full_copy_error_types', 'r')
#     hist = [['d', 0], ['ld', 0], ['s', 0], ['pi', 0], ['n', 0]]
#     lines = full_copy_err_type_f.readlines()
#     for line in lines:
#         if line == 'd\n':
#             hist[0][1] += 1
#         if line == 'ld\n':
#             hist[1][1] += 1
#         if line == 's\n':
#             hist[2][1] += 1
#         if line == 'pi\n':
#             hist[3][1] += 1
#         if line == 'n\n':
#             hist[4][1] += 1
#     full_copy_err_type_f.close()
#
#     full_copy_err_type_ana_f = open('error_types_analysis', 'w')
#     full_copy_err_type_ana_f.write('d appearance rate: ' + str(hist[0][1] / (1000 * len(example_strand))) + '\n')
#     full_copy_err_type_ana_f.write('ld appearance rate: ' + str(hist[1][1] / (1000 * len(example_strand))) + '\n')
#     full_copy_err_type_ana_f.write('s appearance rate: ' + str(hist[2][1] / (1000 * len(example_strand))) + '\n')
#     full_copy_err_type_ana_f.write('pi appearance rate: ' + str(hist[3][1] / (1000 * len(example_strand))) + '\n')
#     full_copy_err_type_ana_f.write('n appearance rate: ' + str(hist[4][1] / (1000 * len(example_strand))) + '\n')
#     full_copy_err_type_ana_f.close()
#
#     error_loc_f = open('error_locations', 'r')
#     hist = [['1/3 len', 0], ['rest', 0]]
#     lines = error_loc_f.readlines()
#     for line in lines:
#         val = line.rstrip()
#         val = int(val)
#         if val <= len(example_strand)/3:
#             hist[0][1] += 1
#         else:
#             hist[1][1] += 1
#
#     error_loc_f.close()
#
#     err_loc_ana_f = open('error_locations_analysis', 'w')
#     err_loc_ana_f.write('1/3 appearance rate: ' + str(hist[0][1] / 1000) + '\n')
#     err_loc_ana_f.write('after 1/3 appearance rate: ' + str(hist[1][1] / 1000) + '\n')
#
#     err_loc_ana_f.close()
#
#     deletion test: - OBSOLETE!
#
#     simulator.err_type = 'd'
#     new_strand = simulator.inject_deletion('d')
#     deletion_f = open('deletion', 'w')
#     deletion_f.write('single base:\n')
#     deletion_f.write('original:\n' + example_strand + '\n')
#     deletion_f.write('modified:\n' + new_strand + '\n')
#
#     deletion_f.write('\n')
#
#     simulator.err_type = 'ld'
#     new_strand = simulator.inject_deletion('ld')
#     deletion_f.write('multiple base:\n')
#     deletion_f.write('original:\n' + example_strand + '\n')
#     deletion_f.write('modified:\n' + new_strand + '\n')
#
#     deletion_f.close()


