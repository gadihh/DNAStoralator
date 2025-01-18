'''
Reads Length Analyzer
Preforming this analysis will result with an additional section in the generated report.
The additional section will include a grphical representation of the distribution
of the length of the reads.
'''

import matplotlib.pyplot as plt
import numpy as np

from src.analyzers.analyzer import Analyzer
from src.utils.content import Content


class ReadsLengthAnalyzer(Analyzer):
    name = "Read length analyzer"
    display_rank = 50
    def analyze(self, library_reads, library_design):
        if not library_reads.did_matching():
            print("You can not preform variant distribution analysis on a library that did not do matching.")
            return
        des_len=library_design.get_design_longest_sequence_len()
        sequences = library_reads.get_original_sequences()
        c_data = []
        seq = library_reads.get_matched_reads()
        for read in seq:
            c_data.append(round(len(str(read.get_attribute('sequence')))))
        print(np.where(c_data is None))
        return self.generate_content(c_data, des_len)

    @staticmethod
    def generate_content(data, des_len):
        plt.xlabel("Length", fontsize=20)
        plt.ylabel("Fraction of reads", fontsize=20)
        plt.yscale('log')

        his, bini = np.histogram(data)

        bins = np.arange(min(data), max(data) + 1.5) -0.5

        # then you plot away
        fig, ax = plt.subplots()
        _ = ax.hist(data, bins, rwidth=0.6, density=True)
        ax.set_xticks(bins + 0.5)
        plt.yscale('log')
        plt.xlabel("Length", fontsize=20)
        plt.ylabel("Fraction of reads", fontsize=20)
        plt.xticks((min(data), int((max(data)+min(data)) / 2), max(data)), [min(data), int((max(data)+min(data)) / 2),
                                                                                      max(data)], rotation=0, fontsize=16)

        plt.yticks(fontsize=16)
        plt.grid()
        imgname="temp/reads_length_histogram.png"

        plt.savefig(imgname, dpi=300, bbox_inches = "tight")
        plt.clf()

        return [Content(Content.Type.TEXT, "Reads length distribution"),
                Content(Content.Type.IMAGE, imgname)]
