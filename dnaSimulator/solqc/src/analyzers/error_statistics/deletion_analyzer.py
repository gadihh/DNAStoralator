from __future__ import print_function
from progress.bar import Bar

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import edlib
import re

from dnaSimulator.solqc.src.analyzers.analyzer import Analyzer
from dnaSimulator.solqc.src.utils.content import Content
from dnaSimulator.solqc.src.utils.biology import DELETION, INSERTION, MISMATCH


class DeletionAnalyzer (Analyzer):
    name = "Deletion Analyzer"
    requires_alignment = True

    # Will store the number of deletions for each base.
    base_deletion = {}

    # Stores the number of deletion in each position
    position_deletion = []

    def __init__(self):
        super(DeletionAnalyzer, self).__init__()
        self.base_deletion = {
            'A' : 0,
            'G' : 0,
            'C' : 0,
            'T' : 0
        }

    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance():
            print("You can prefrom deletion analyzing on reads that did not go through edit distance")
            return None

        # Initalize the position deletion to the size of the longest read in the design
        # TODO each position should be divided by the amount of times it appears (position 1 appears much more then 150)
        self.position_deletion = np.zeros(library_design.get_design_longest_sequence_len())

        # Get progress bar object
        bar = self.get_progress_bar(library_reads.get_matched_reads_count())

        for read in library_reads.iterate_matched_reads():
            # TODO find a better name.
            query_target_path = read.get_cigar_path()
            read_count = read.get_row_count()

            # Update Deletion counts
            self.update_deletion(query_target_path, read(), read_count)

            # # Update progress
            bar.next()

        bar.finish()
        return self.generate_content()

    def generate_content(self):
        plt.xlabel("position")
        plt.ylabel("Number of deletions")
        sns.regplot(np.arange(len(self.position_deletion)), self.position_deletion)
        file_name = "temp/deletion_count.png"
        plt.savefig(file_name)
        plt.clf()

        return [Content(Content.Type.TEXT, "Deletion Count"),
                Content(Content.Type.IMAGE, file_name)]

    def update_deletion(self, path, variant, read_count):
        deletion_indices = locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            if i < len(self.position_deletion):
                self.position_deletion[i] += read_count

    @staticmethod
    def locate_deletion_locations(path):
        offset = 0
        indices = []
        for index, letter in enumerate(str(path)):
            if letter == DELETION:
                indices.append(index - offset)
            if letter == INSERTION:
                offset += 1

        return indices

    def conclude(self):
        sns.regplot(np.arange(len(self.position_deletion)), self.position_deletion)
        # plt.xlabel("position")
        # plt.ylabel("Number of deletions")
        plt.show()

    def __str__(self):
        return "Deletion"


def locate_deletion_locations(path):
    offset = 0
    indices = []

    for index, letter in enumerate(path):
        if letter == DELETION:
            indices.append(index - offset)
        if letter == INSERTION:
            offset += 1

    return indices


if __name__ == "__main__":
    query = "CCAATATCCTTAGCTGATCACCCATATGCCCACTAAGTACCCGGAAGAGCGAATGGCGTAGTGCCATGTTCCTTATTGCCCGGGCCACAGTTGACATTAGATATGG"
    target = "CCAATATCCTTAGCTGATCACCCATATGCCCACTAAGTACCCGGAAGAGCGAATGGCGTAGTGCCGCAGCCACACACAGCAGGAGCAGAATTCGTGTCAACATGGGCATGTTCCTTATTGCCCGGGCCACAGTTGACATTAGATATGG"
    align = edlib.align(query,
                        target,
                        task="path")
    print("Query length = {}\t Target length = {}".format(len(query), len(target)))
