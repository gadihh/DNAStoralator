import matplotlib.pyplot as plt
import numpy as np

from src.analyzers.analyzer import Analyzer
from src.utils.content import Content


class GeneralAnalyzer(Analyzer):
    display_rank = 10
    name = "General analyzer"

    def analyze(self, library_reads, library_design):
        if not library_reads.did_matching():
            print("You can not preform variant distribution analysis on a library that did not do matching.")
            return

        matched_sequences_num = library_reads.get_matched_reads_count()
        unmatched_sequences_num = library_reads.get_unmatched_reads_count()
        original_sequences_num = matched_sequences_num + unmatched_sequences_num

        content_data = {
            "original_sequences_num" : original_sequences_num,
            "matched_sequences_num" : matched_sequences_num,
            "unmatched_sequences_num" : unmatched_sequences_num
        }
        return self.generate_content(content_data)

    @staticmethod
    def generate_content(data):
        original_sequences_num = data['original_sequences_num']
        matched_sequences_num = data['matched_sequences_num']
        unmatched_sequences_num = data['unmatched_sequences_num']
        past_preprocess = matched_sequences_num + unmatched_sequences_num

        text_1 = "Started with {} reads".format(original_sequences_num)
        text_2 = "Number of reads after pre-processing = {}".format(matched_sequences_num + unmatched_sequences_num)
        text_3 = "Fraction of sequences who passed pre-processing = {0:.2f}".format(past_preprocess / original_sequences_num)

        text_4 = "Number of matched reads = {}".format(matched_sequences_num)
        text_5 = "Fraction of sequences who got matched from pre-processed= {0:.2f}".format(matched_sequences_num / past_preprocess)
        text_6 = "Fraction of sequences who got matched from all reads= {0:.2f}".format(matched_sequences_num / original_sequences_num)

        return [Content(Content.Type.TEXT, text_1),
                Content(Content.Type.TEXT, text_2),
                Content(Content.Type.TEXT, text_3),
                Content(Content.Type.TEXT, text_4),
                Content(Content.Type.TEXT, text_5),
                Content(Content.Type.TEXT, text_6)]