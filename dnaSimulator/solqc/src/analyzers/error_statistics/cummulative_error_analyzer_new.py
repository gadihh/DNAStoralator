from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from src.analyzers.error_statistics.deletion_analyzer import DeletionAnalyzer
from src.utils.content import Content
from matplotlib.pyplot import axes as ax


DELETION = "1"
INSERTION = "2"
MISMATCH = "3"

CIGAR_DICT = {
    "=": "0",
    "I": INSERTION,
    "D": DELETION,
    "X": MISMATCH
}


def plot_errornumcount(x_labels, y_values, graph_title='', x_title='', y_title='', x_title_fontsize=14, y_title_fontsize=14,
                       x_labels_fontsize=14, multigraph_labels=[''], bar_width=0.15, plot_format=['g-o', 'b', 'r', 'y'],
                       indexing_jump=1, maxxtick=45):
    fig = plt.figure(figsize=(8, 4), linewidth=3)
    # index = np.arange(len(x_labels))
    index= np.arange(maxxtick)
    #print(x_labels[1])
    #plt.yscale('log')


    opacity = 0.8

    for i, values in enumerate(y_values):
        rects1 = plt.plot(x_labels[0:maxxtick], y_values[i][0:maxxtick], plot_format[i], label=multigraph_labels[i])
    plt.xticks(list(reversed(list(reversed(index))[0::indexing_jump]))[0:maxxtick],
               list(reversed(list(reversed(x_labels[0:maxxtick]))[0::indexing_jump])),
               rotation=0, fontsize=36)  # index + bar_width?

    plt.yticks(fontsize=36)
    #plt.gca().set_yticklabels(['{:.0f}%'.format(x * 100) for x in plt.gca().get_yticks()])
    plt.tight_layout()
    plt.grid()
    plt.xlabel("Number of edit errors", fontsize=36)
    plt.ylabel("Fraction of reads", fontsize=36)
    #locs, labels = plt.xticks()
    #plt.xticks(locs[::5], labels[::5])
    # Save scatter as image
    image_one_name = 'temp/cummulativeDistribution.png'
    plt.yscale('log')
    plt.savefig(image_one_name, dpi=300, bbox_inches = "tight")
    #plt.show()

    plt.clf()

    return image_one_name

    # plt.show()


def plot_errornumcount_with_long(x_labels, y_values, graph_title='', x_title='', y_title='', x_title_fontsize=14,
                                 y_title_fontsize=14,
                                 x_labels_fontsize=14, multigraph_labels=[''], bar_width=0.15,
                                 plot_format=['g-o', 'b', 'r', 'y'],
                                 indexing_jump=1, maxxtick=45):
    fig = plt.figure(figsize=(8, 4), linewidth=3)
    # index = np.arange(len(x_labels))
    index = np.arange(maxxtick)
    print(x_labels[1])
    # plt.yscale('log')

    opacity = 0.8

    for i, values in enumerate(y_values):
        rects1 = plt.plot(x_labels[0:maxxtick], y_values[i][0:maxxtick], plot_format[i], label=multigraph_labels[i],
                          linewidth=2.0)

    plt.xticks(list(reversed(list(reversed(index))[0::indexing_jump]))[0:maxxtick],
               list(reversed(list(reversed(x_labels[0:maxxtick]))[0::indexing_jump])),
               rotation=0, fontsize=36)  # index + bar_width?
    # plt.legend(loc=3, fontsize=16)
    # plt.ylim(0, 0.12)
    plt.yticks(fontsize=36)
    plt.gca().set_yticklabels(['{:.0f}%'.format(x * 100) for x in plt.gca().get_yticks()])
    plt.tight_layout()
    plt.grid()
    plt.xlabel("Number of edit errors", fontsize=36)
    plt.ylabel("Fraction of reads", fontsize=36)
    #locs, labels = plt.xticks()
    #plt.xticks(locs[::5], labels[::5])
    # Save scatter as image
    image_one_name = 'temp/inputWIthXorLessError_withLong.png'.format('A')  # TODO: change the pic name :)
    plt.savefig(image_one_name, dpi=300, bbox_inches = "tight")
    #plt.show()
    plt.clf()

    return image_one_name

    # plt.show()



class NewCumulativeErrorAnalyzer(DeletionAnalyzer):
    # Will store the number of deletions for each base.
    base_deletion = {}

    def __init__(self):
        super(DeletionAnalyzer, self).__init__()
        self.read_errors_count = []
        self.read_errors_count_with_long = []
        self.longest_read_length = 0

        for key in self.base_deletion:
            self.base_deletion[key] = np.zeros(self.longest_read_length)

    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance:
            print("You can prefrom deletion analyzing on reads that did not go through edit distance")
            return None
        self.longest_read_length = library_design.get_design_longest_sequence_len()

        # Get all the reads we managed to match.
        matched_reads = library_reads.get_matched_reads()
        bar = self.get_progress_bar(len(matched_reads))
        max_error_counter=0
        max_error_counter_long=0
        for read in matched_reads:
            query_target_path=read.get_cigar_path()
            current_counter=self.get_count_error(query_target_path)
            #if current_counter >= 1000:
            #    print(read.get_attribute('sequence'))
            #    print(read.get_variant_id())
            #    print(query_target_path)
            #    exit(0)
            current_counter_long=self.get_count_error_with_long(query_target_path)
            if(current_counter>max_error_counter):
                max_error_counter=current_counter
            if(current_counter_long>max_error_counter_long):
                max_error_counter_long=current_counter_long
        for i in range(0, max_error_counter):
            self.read_errors_count.append(0)
        for i in range(0, max_error_counter_long):
            self.read_errors_count_with_long.append(0)
        print("symbol*****************")
        print(max_error_counter)
        print(max_error_counter_long)

        for read in matched_reads:
            # TODO find a better name.
            query_target_path = read.get_cigar_path()
            variant = library_design.get_variant_by_read(read)

            # count the deletion number of the read
            self.count_error_num(query_target_path, max_error_counter, read.get_row_count())
            self.count_error_num_with_long(query_target_path, max_error_counter_long, read.get_row_count())
            # Update Deletion counts

            # # Update progress
            bar.next()

        bar.finish()

        return self.generate_content2(self.read_errors_count, library_reads.get_matched_reads_count())

    def __str__(self):
        return 'Per base deletion analyzer'

    def generate_content2(self, read_errors_count, len_matched_reads):
        content_array = []
        colors = {
            "A": "#0099ff",
            "C": "#ff9900",
            "G": "#7fff99",
            "T": "#e60042"
        }

        headline = Content(Content.Type.TEXT, "% of reads with X or less errors")
        content_array.append(headline)
        for idx, item in enumerate(read_errors_count):
            read_errors_count[idx] = (item / len_matched_reads)
        dists = np.zeros((1, len(read_errors_count))).astype(float)
        dists[0] = read_errors_count
        multigraph_labels = ['E', 'C', 'T', 'G']

        image_one_name = plot_errornumcount(x_labels=np.array(range(0, max(11, int(len(read_errors_count)/8)))), y_values=dists,
                                            x_title='Index',
                                            y_title='Probability', graph_title='Error Distribution Per Index',
                                            multigraph_labels=multigraph_labels, indexing_jump=1, x_labels_fontsize=20,
                                            maxxtick=max(11, int(len(read_errors_count)/8)))
        # Create and add content
        image_one = Content(Content.Type.IMAGE, image_one_name)

        content_array.append(image_one)

        #headline = Content(Content.Type.TEXT, "% of reads with X or less errors")
        content_array.append(headline)
        for idx, item in enumerate(self.read_errors_count_with_long):
            self.read_errors_count_with_long[idx] = (item / len_matched_reads)
        dists = np.zeros((1, len(self.read_errors_count_with_long))).astype(float)
        dists[0] = self.read_errors_count_with_long
        image_one_name = plot_errornumcount_with_long(
            x_labels=np.array(range(0, max(11, int(len(self.read_errors_count_with_long)/8)))), y_values=dists, x_title='Index',
            y_title='Probability', graph_title='Error Distribution Per Index',
            multigraph_labels=multigraph_labels, indexing_jump=1, x_labels_fontsize=20, maxxtick=max(11, int(len(self.read_errors_count_with_long)/8)))
        image_one = Content(Content.Type.IMAGE, image_one_name)

        content_array.append(image_one)

        return content_array

    def update_deletion(self, path, variant):
        deletion_indices = self.locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            self.base_deletion[variant[i]][i] += 1

    def count_error_num(self, path, max_error_counter, count):
        counter = 0
        #if float('-inf') < float(path) < float('inf'):
        if True:
            for index, letter in enumerate(path):
                if letter not in {'0', '1', '2', '3'}:
                    print(letter)
                    return counter
            for index, letter in enumerate(path):
                if letter not in {'0', '1', '2', '3'}:
                    print(letter)
                    break
                if letter != '0':
                    if letter == '1':
                        if index-1 >= 0 and str(path)[index-1] != '1':
                            counter+=count
                        #while index+1 < len(str(path)) and str(path)[index+1] == '1':
                        #    index+=1
                        #counter+=count
                    else:
                        counter+=count
                        #while index + 1 < len(str(path)) and str(path)[index + 1] == '1':
                        #    index += 1
                        #counter += count

        for j in range(counter, max_error_counter):
            # print(max_error_counter-j)
            # self.read_errors_count[max_error_counter-j-1]=self.read_errors_count[max_error_counter-j-1]+1;
            self.read_errors_count[j] = self.read_errors_count[j] + count;

    def get_count_error(self, path):
        counter = 0
        #if float('-inf') < float(path) < float('inf'):
        if True:
            for index, letter in enumerate(path):
                if letter not in {'0', '1', '2', '3'}:
                    break
                if letter != '0':
                    if letter == '1':
                        if index-1 >= 0 and str(path)[index-1] != '1':
                            counter+=1
                        #while index+1 < len(str(path)) and str(path)[index+1] == '1':
                        #    index+=1
                        #counter+=1
                    else:
                        counter+=1
                        #while index + 1 < len(str(path)) and str(path)[index + 1] == '1':
                        #    index += 1
                        #counter+=1

        return counter


    def get_count_error_with_long(self, path):
        counter = 0
        #if float('-inf') < float(path) < float('inf'):
        if True:
            for index, letter in enumerate(path):
                if letter not in {'0', '1', '2', '3'}:
                    print(letter)
                    break
                if letter != '0':
                    counter += 1
            return counter
        else:
            return counter
        endif


    def count_error_num_with_long(self, path, max_error_counter, read_count):
        counter = 0
        #if float('-inf') < float(path) < float('inf'):
        if True:
            for index, letter in enumerate(path):
                if letter not in {'0', '1', '2', '3'}:
                    print(letter)
                    break
                if letter != '0':
                    counter+=read_count
        for j in range(counter, max_error_counter):
            #print(max_error_counter-j)
            #self.read_errors_count[max_error_counter-j-1]=self.read_errors_count[max_error_counter-j-1]+1;
            self.read_errors_count_with_long[j]=self.read_errors_count_with_long[j]+read_count;

    def get_position_base_count(self, library_reads, library_design):
        """
        Compute how many bases should ideally be in each position if all the reads where perfect.
        :param library_reads:
        :param library_design:
        :return: A dictionary with all the letters and for each letter an array where each index
        corresponds to a position and each value indicates how many of the base should be in that position.
        """

        # Starting each array with ones,
        letter_position = {
            'A': np.ones(self.longest_read_length),
            'G': np.ones(self.longest_read_length),
            'T': np.ones(self.longest_read_length),
            'C': np.ones(self.longest_read_length)
        }

        df_copy = library_reads.get_df_copy()
        group_by_variant = df_copy.groupby('variant_id')

        for variant_id, group in group_by_variant:
            if variant_id == -1:
                continue

            num_of_reads = len(group)
            variant = library_design.get_variant_by_id(variant_id)
            # Converting variant sequence to a numpy string array for vectorization purposes.
            v_sequence = np.asarray(list(variant()))
            for letter in letter_position:
                letter_in_sequence = (v_sequence == letter)
                letter_in_sequence = letter_in_sequence * num_of_reads

                # Some variants are shorter which will cause a bug. This fixes the bug. yeah!
                sequence_length = len(v_sequence)
                letter_position[letter][:sequence_length] += letter_in_sequence

        return letter_position
