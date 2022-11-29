import argparse
import threading
from leven import levenshtein

from simulator import *
import time


class Clustering:
    def __init__(self, chosen_index, chosen_technology):
        self.index = chosen_index
        self.technology = chosen_technology
        self.letters = ['A', 'C', 'G', 'T']
        self.dict_indices = {}
        self.dict_strings = {}
        self.evyat_dict_strings = {}
        self.thrown_strands_dict = {}
        if platform.system() == "Linux":
            self.shuffled_file = '/home_nfs/sgamer8/DNAindex' + str(
                self.index) + '/files/' + self.technology + '/' + 'errors_shuffled.txt'
        elif platform.system() == "Windows":
            self.shuffled_file = 'files/' + self.technology + '/' + 'errors_shuffled.txt'
        self.number_of_strands_in_wrong_cluster = 0
        self.number_of_false_negative = 0
        self.total_amount_of_thrown = 0

    # creates dictionaries with all the possible keys of length 5, and an empty list of values
    def create_all_keys(self):
        if self.index == 9:
            for first in self.letters:
                for second in self.letters:
                    for third in self.letters:
                        for fourth in self.letters:
                            for fifth in self.letters:
                                for sixth in self.letters:
                                    for seventh in self.letters:
                                        for eight in self.letters:
                                            for ninth in self.letters:
                                                possible_index = first + second + third + fourth + fifth + sixth + seventh + eight + ninth
                                                self.dict_indices[possible_index] = []
                                                self.dict_strings[possible_index] = []
                                                self.evyat_dict_strings[possible_index] = []
                                                self.thrown_strands_dict[possible_index] = []
        elif self.index == 5:
            for first in self.letters:
                for second in self.letters:
                    for third in self.letters:
                        for fourth in self.letters:
                            for fifth in self.letters:
                                possible_index = first + second + third + fourth + fifth
                                self.dict_indices[possible_index] = []
                                self.dict_strings[possible_index] = []
                                self.evyat_dict_strings[possible_index] = []
                                self.thrown_strands_dict[possible_index] = []
        elif self.index == 4:
            for first in self.letters:
                for second in self.letters:
                    for third in self.letters:
                        for fourth in self.letters:
                            possible_index = first + second + third + fourth
                            self.dict_indices[possible_index] = []
                            self.dict_strings[possible_index] = []
                            self.evyat_dict_strings[possible_index] = []
                            self.thrown_strands_dict[possible_index] = []
        print('Finished creating the dictionary with the possible keys')

    # fills the dictionaries with values from the shuffled file (one dict with indices and one with the actual strings)
    # also outputs the dictionaries into files. the indices into the allclustersdict file and the strings into
    # allclustersstringdict file
    def fill_dict_from_shuffled(self):
        with open(self.shuffled_file, 'r') as evyat_shuffled:
            for counter, line in enumerate(evyat_shuffled):
                # print(counter, line[:5])
                self.dict_indices[line[:self.index]].append(counter)
                self.dict_strings[line[:self.index]].append(line.rstrip('\n'))
        if platform.system() == "Linux":
            outputdict_path = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/cluster_output/' + self.technology + '/' + str(
                self.index) + '_allclustersdict'
            outputstringsdict_path = '/home_nfs/sgamer8/DNAindex' + str(
                self.index) + '/cluster_output/' + self.technology + '/' + str(self.index) + '_allclustersstringdict'
        elif platform.system() == "Windows":
            outputdict_path = 'cluster_output/' + self.technology + '/' + str(self.index) + '_allclustersdict'
            outputstringsdict_path = 'cluster_output/' + self.technology + '/' + str(self.index) + '_allclustersstringdict'
        with open(outputdict_path, 'w') as outputdict, open(outputstringsdict_path, 'w') as outputstringsdict:
            for key1, key2 in zip(self.dict_indices, self.dict_strings):
                print(key1, end=" ", file=outputdict)
                print(self.dict_indices[key1], file=outputdict)
                print(key2, end=" ", file=outputstringsdict)
                print(self.dict_strings[key2], file=outputstringsdict)
        print('Finished filling the dictionary with keys from the shuffled file')

    # cleans the file of random strings taken from the website
    # taken from this website: https://users-birc.au.dk/palle/php/fabox/random_sequence_generator.php
    def clean_random_file(self):
        with open('files/256_seqs_196_196_bp.fasta', 'r') as inputfile:
            with open('files/clean_256_seqs_196_196_bp.fasta', 'w') as outputfile:
                for line in inputfile:
                    if 'random sequence' in line.strip():
                        continue
                    elif line.strip() == '':
                        continue
                    else:
                        outputfile.write(line.upper())

    # sticks the index for each strand in the beginning of the random clean strands file by order
    def unite_letters_to_clean_file(self, index_size):
        with open('files/clean_256_seqs_196_196_bp.fasta', 'r') as clean_file, open(
                'files/' + str(index_size) + '_letters.txt',
                'r') as letters, open(
            'files/' + str(index_size) + '_allclustersofindex.txt', 'w', newline='') as outputfile:
            for line1, line2 in zip(letters, clean_file):
                outputfile.write(line1[:index_size] + line2)

    def create_evyat_dict(self):
        in_cluster = 0

        if platform.system() == "Linux":
            evyat_path = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/files/' + self.technology + '/' + 'evyat.txt'
            evyatdict_path = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/files/' + self.technology + '/' + str(
                self.index) + '_evyatdict'
        elif platform.system() == "Windows":
            evyat_path = 'files/' + self.technology + '/' + 'evyat.txt'
            evyatdict_path = 'files/' + self.technology + '/' + str(self.index) + '_evyatdict'

        with open(evyat_path, 'r') as evyatfile:
            for line in evyatfile:
                if line.strip() != '' and '*' not in line.strip() and in_cluster == 0:
                    in_cluster = 1
                    cluster_index = line[:self.index]
                    continue
                if line.strip() == '':
                    in_cluster = 0
                    continue
                if '*' in line.strip():
                    continue
                if line.strip() != '' and '*' not in line.strip() and in_cluster == 1:
                    self.evyat_dict_strings[cluster_index].append(line.rstrip('\n'))

        with open(evyatdict_path, 'w') as outputdict:
            for key in self.evyat_dict_strings:
                print(key, end=" ", file=outputdict)
                print(self.evyat_dict_strings[key], file=outputdict)

    def compare_evyat_with_clustering(self):
        if platform.system() == "Linux":
            strands_in_wrong_cluster = '/home_nfs/sgamer8/DNAindex' + str(
                self.index) + '/cluster_output/' + self.technology + '/' + str(self.index) + '_strands_in_wrong_cluster.txt'
            false_positive = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/cluster_output/' + self.technology + '/' + str(
                self.index) + '_false_positive.txt'
        elif platform.system() == "Windows":
            strands_in_wrong_cluster = 'cluster_output/' + self.technology + '/' + str(
                self.index) + '_strands_in_wrong_cluster.txt'
            false_positive = 'cluster_output/' + self.technology + '/' + str(self.index) + '_false_positive.txt'
        with open(strands_in_wrong_cluster, 'w') as wrong_strands_file, open(false_positive,
                                                                             'w') as false_positive_file:
            for key in self.dict_strings.keys():
                # TODO print the full key in the new evyat file (maybe this one doesn't have the full one, so might need to make a dict between the index and the strand
                # TODO then print *****************************
                real_cluster_for_key = self.evyat_dict_strings[key]
                output_cluster_for_key = self.dict_strings[key]
                for strand in output_cluster_for_key:
                    if strand not in real_cluster_for_key:
                        self.number_of_strands_in_wrong_cluster = self.number_of_strands_in_wrong_cluster + 1
                        print('Strand: ' + strand, file=wrong_strands_file)
                        print('In the wrong cluster: ' + key, file=wrong_strands_file)
                        print('***************************', file=wrong_strands_file)
                    # TODO else: write it to the new evyat
                for strand in self.thrown_strands_dict[key]:
                    if strand in real_cluster_for_key:
                        self.number_of_false_negative = self.number_of_false_negative + 1
                        print('Thrown strand: ' + strand, file=false_positive_file)
                        print('Cluster it should be in: ' + key, file=false_positive_file)
                        print('***************************', file=false_positive_file)
                # TODO print 2 new lines and then go to the next key
            print('Number of strands in wrong cluster: ' + str(self.number_of_strands_in_wrong_cluster))
            print('Number of false negative: ' + str(self.number_of_false_negative))
            print('Total number of thrown strands: ' + str(self.total_amount_of_thrown))
            return [str(self.number_of_strands_in_wrong_cluster), str(self.number_of_false_negative)]

    # creates a file with all the indices of size 5 or 9.
    # default is 9
    def create_indices_only(self, index_size=9):
        if index_size == 9:
            with open('files/9_letters.txt', 'w') as index_file:
                for first in self.letters:
                    for second in self.letters:
                        for third in self.letters:
                            for fourth in self.letters:
                                for fifth in self.letters:
                                    for sixth in self.letters:
                                        for seventh in self.letters:
                                            for eight in self.letters:
                                                for ninth in self.letters:
                                                    possible_index = first + second + third + fourth + fifth + sixth + seventh + eight + ninth
                                                    index_file.write(possible_index + '\n')
        elif index_size == 5:
            with open('files/5_letters.txt', 'w') as index_file:
                for first in self.letters:
                    for second in self.letters:
                        for third in self.letters:
                            for fourth in self.letters:
                                for fifth in self.letters:
                                    possible_index = first + second + third + fourth + fifth
                                    index_file.write(possible_index + '\n')
        elif index_size == 4:
            with open('files/4_letters.txt', 'w') as index_file:
                for first in self.letters:
                    for second in self.letters:
                        for third in self.letters:
                            for fourth in self.letters:
                                possible_index = first + second + third + fourth
                                index_file.write(possible_index + '\n')

    def full_edit_distance_threads(self):
        if platform.system() == "Linux":
            edit_dist_avg_f = '/home_nfs/sgamer8/DNAindex9/output/9_edit_dist_averages.txt'
        elif platform.system() == "Windows":
            edit_dist_avg_f = 'output/5_edit_dist_averages.txt'
        with open(edit_dist_avg_f, 'w') as avg_file:
            threads = []
            for key in list(self.dict_strings.keys()):
                t = threading.Thread(target=self.full_edit_distance_per_cluster, args=(key,))
                t.start()
                threads.append(t)
                print(key)
                # print(key + ' - ' + total_avg, file=avg_file)
            for t in threads:
                t.join()

    def full_edit_distance(self, report_func):
        if platform.system() == "Linux":
            edit_dist_avg_f = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/cluster_output/' + self.technology + '/' + str(
                self.index) + '_edit_dist_averages.txt'
        elif platform.system() == "Windows":
            edit_dist_avg_f = 'cluster_output/' + self.technology + '/' + str(self.index) + '_edit_dist_averages.txt'

        num_values = pow(4, self.index)
        i = 0

        with open(edit_dist_avg_f, 'w') as avg_file:
            for key in list(self.dict_strings.keys()):
                report_func(num_values, i)
                total_avg = self.full_edit_distance_per_cluster(key)
                i = i + 1
                print(key)
                print(key + ' - ' + total_avg, file=avg_file)

    def full_edit_distance_per_cluster(self, givenkey):

        levenstein_dict = []
        averages_dict = []
        averages_list = []
        output_cluster_for_key = self.dict_strings[givenkey]
        # for strand in output_cluster_for_key:
        #     print (strand)
        # print('--------------------------------------------')
        for strand in output_cluster_for_key:
            levenstein_dict.append({strand: []})
            averages_dict.append({strand: []})
            # levenstein_dict[strand] = []
            # averages_dict[strand] = []
        for counter1, strand in enumerate(output_cluster_for_key):
            for counter2, strand2 in enumerate(output_cluster_for_key):
                edit_distance = levenshtein(strand, strand2)
                levenstein_dict[counter1][strand].append(edit_distance)
        for counter, liststrand in enumerate(levenstein_dict):
            averages_dict[counter][list(liststrand.keys())[0]].append(Average(levenstein_dict[counter][list(liststrand.keys())[0]]))
            averages_list.append(averages_dict[counter][list(liststrand.keys())[0]])
            # print('Avg: ' + str(averages_dict[key]) + '---' + key + ': ' + str(levenstein_dict[key]))
        final_avg_list = []
        for avg in averages_list:
            final_avg_list.append(avg[0])
        total_average = Average(final_avg_list) + 10
        # print('Total average: ' + str(total_average))
        for counter, dict in enumerate(averages_dict):
            if averages_dict[counter][list(dict.keys())[0]][0] > total_average:
                self.thrown_strands_dict[givenkey].append(list(dict.keys())[0])
                output_cluster_for_key.remove(list(dict.keys())[0])
                self.total_amount_of_thrown = self.total_amount_of_thrown + 1
            # print(averages_dict[key])
        self.dict_strings[givenkey] = output_cluster_for_key
        # outputcluster = self.dict_strings['AAAAA']
        # for strand in outputcluster:
        #     print (strand)
        return str(total_average)


def Average(lst):
    return sum(lst) / len(lst)


def parse_cmd_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--technology", type=str, default='none', help="selected technology for simulator")
    parser.add_argument("-i", "--index", type=int, default=5, help="index size to be used")
    parser.add_argument("-s", "--simulator_status", type=str, default='off', help="simulator on or off")
    args = parser.parse_args()
    return args.technology, args.index, args.simulator_status

