'''
Frequency Analyzer
Preforming this analysis will result with an outputted frequecny_count.csv file.
The frequency count will contain 2 columns : variant id, number of reads.
variant_id : The row number of the variant in the design file.
number_of_reads : The number of reads found for the corresponding variant.
'''

from src.analyzers.analyzer import Analyzer
import src.config as config


class FrequencyAnalyzer(Analyzer):
    name = "Frequency Analyzer"

    def analyze(self, library_reads, library_design):
        if not library_reads.did_matching():
            print("You can not preform variant frequency analysis on a library that did not do matching.")
            return

        reads_df = library_reads.get_df_copy()

        def get_variant_sequence(row, ld):
            return ld.get_variant_sequence(row.name)

        frequency_df = reads_df.groupby(['variant_id'])['count'].sum().to_frame('count')
        frequency_df['variant_sequence'] = frequency_df.apply(get_variant_sequence, ld=library_design, axis=1)

        save_path = config.get_deliverable_path()
        frequency_df.to_csv('{}/frequency_count.csv'.format(save_path))
