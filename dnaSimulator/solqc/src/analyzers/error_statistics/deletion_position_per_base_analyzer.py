"""
Deletion per position per base analyzer:
Preforming this analysis will result with an additional section in the generated report.
This will display the amount of deletion in each position.
1 plot to show all deletion, and 4 plots for each letter.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from src.analyzers.error_statistics.deletion_analyzer import DeletionAnalyzer
from src.utils.content import Content


class DeletionPositionPerBaseAnalyzer (DeletionAnalyzer):
    name = "Deletion position per base analyzer"
    requires_alignment = True

    # Will store the number of deletions for each base.
    base_deletion = {}

    def __init__(self):
        super(DeletionAnalyzer, self).__init__()
        self.base_deletion = {
            'A': 0,
            'G': 0,
            'C': 0,
            'T': 0,
        }

        self.longest_read_length = 0

    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance():
            print("You can not perform deletion analyzing on reads that did not go through edit distance")
            return None

        self.longest_read_length = library_design.get_design_longest_sequence_len()

        for key in self.base_deletion:
            self.base_deletion[key] = np.zeros(self.longest_read_length)

        # Get all the reads we managed to match.
        matched_reads = library_reads.get_matched_reads()

        bar = self.get_progress_bar(len(matched_reads))
        for read in matched_reads:
            # TODO find a better name.
            query_target_path = read.get_cigar_path()
            variant = library_design.get_variant_by_read(read)

            read_count = read.get_row_count()
            # Update Deletion counts
            self.update_deletion(query_target_path, variant(), read_count)

            # # Update progress
            bar.next()

        bar.finish()

        # Dividing each position by the amount of expected base in that position..
        letter_position = self.get_position_base_count(library_reads, library_design)
        letter_position_reads = self.get_base_per_position_by_reads(library_reads, library_design)

        for key, value in self.base_deletion.items():
            self.base_deletion[key] /= letter_position[key]

        return self.generate_content(letter_position, letter_position_reads)

    def __str__(self):
        return 'Per base deletion analyzer'

    def generate_content(self, letter_position, letter_position_reads):
        content_array = []
        colors = {
            "A": "#0099ff",
            "C": "#ff9900",
            "G": "#7fff99",
            "T": "#e60042"
        }

        headline = Content(Content.Type.TEXT, "Base Deletion Ratio Per Position")
        content_array.append(headline)

        for letter, deletion_position_rate in self.base_deletion.items():
            # Create scatter plot
            plt.xlabel("{}".format(letter))
            plt.ylabel("Deletion rate")
            sns.regplot(np.arange(len(deletion_position_rate)), deletion_position_rate, color=colors[letter], fit_reg=False)

            # Save scatter as image
            image_one_name = 'temp/a_d_{}.png'.format(letter)
            plt.savefig(image_one_name)
            plt.clf()

            # Create and add content
            image_one = Content(Content.Type.IMAGE, image_one_name)
            content_array.append(image_one)

        # Generate stacked bar graph for base position distribution.
        df = pd.DataFrame(letter_position)
        df = df[['A', 'G', 'T', 'C']] # data frame includes a column for the position we want to get rid of.
        n_position = len(df)

        image_name = "temp/stacked_letters.png"

        df.plot(kind='bar', stacked=True, legend='reverse')
        plt.xlabel("Position")
        plt.ylabel("Count (Stacked on top each other)")
        # Setting x-ticks otherwise the labels over lap each other and look terrible.
        plt.xticks((0, int(n_position / 2), n_position), [0, int(n_position / 2), n_position], rotation=0)
        plt.title("Base distribution by design")
        #plt.legend()
        plt.xlim(0,n_position+5)

        plt.savefig(image_name)
        plt.clf()

        content_array.append(Content(Content.Type.IMAGE, image_name))

        # Generate stacked bar graph for base position distribution by READS
        df = pd.DataFrame(letter_position_reads)
        image_name = "stacked_letters_by_reads.png"
        df.plot(kind='bar', stacked=True, legend='reverse').get_figure()
        plt.xticks((0, int(n_position / 2), n_position), [0, int(n_position / 2), n_position], rotation=0)
        plt.title("Base distribution by reads")
        plt.xlim(0,n_position+5)
        plt.savefig(image_name)

        plt.clf()
        content_array.append(Content(Content.Type.IMAGE, image_name))

        return content_array

    def update_deletion(self, path, variant, read_count):
        deletion_indices = self.locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            self.base_deletion[variant[i]][i] += read_count


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
            'A' : np.ones(self.longest_read_length),
            'G' : np.ones(self.longest_read_length),
            'T' : np.ones(self.longest_read_length),
            'C' : np.ones(self.longest_read_length)
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

    def get_base_per_position_by_reads(self, library_reads, library_design):
        """
        Compute how many bases appeared in each position for the reads.
        :param library_reads:
        :param library_design:
        :return: A dictionary with all the letters and for each letter an array where each index
        corresponds to a position and each value indicates how many of the base were found in the read in that position.
        """

        # Starting each array with ones,
        letter_position = {
            'A' : np.ones(self.longest_read_length+20),
            'G' : np.ones(self.longest_read_length+20),
            'T' : np.ones(self.longest_read_length+20),
            'C' : np.ones(self.longest_read_length+20)
        }

        #df_copy = library_reads.get_df_copy()
        #group_by_variant = df_copy.groupby('variant_id')

        for sequence in library_reads:
            if sequence == -1:
                continue

            # Converting sequence to a numpy string array for vectorization purposes.
            v_sequence = np.asarray(list(sequence()))

            for idx, let in enumerate(v_sequence):
                if idx < len(letter_position[let]):
                  letter_position[let][idx]+=1

        return letter_position
