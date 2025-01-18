'''
Preforming this analysis will result with an outputted reads_matching.csv file.
This csv will contain 3 columns : Read, Count, variant_id(indicating the line of matched variant in the design file).
'''

from src.analyzers.analyzer import Analyzer
import src.config as config


class MatchingAnalyzer(Analyzer):
    name="MatchingAnalyzer"
    # display_rank = 0
    # requires_alignment = False

    def analyze(self, library_reads, library_design):
        reads_matching_df = library_reads.get_df_copy()

        save_path = config.get_deliverable_path()
        reads_matching_df.to_csv('{}/reads_matching.csv'.format(save_path), index=False)
        print("Saved frequency analyzer to : {}".format(save_path))