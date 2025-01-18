# coding=utf-8
"""
This class represents the Oligo library.
We'll keep a reference for all the sequences and will enable easy quereing.
"""
import edlib
import pandas as pd
import json
import swifter

from dnaSimulator.solqc.src.library_reads.read import Read
from dnaSimulator.solqc.src.utils import biology

SEQUENCES_COLUMN = "sequence"
COUNT_COLUMN = "count"
ID_COLUMN = "variant_id"



class LibraryReads(object):
    config = False
    original_sequences = []

    def __init__(self, data, library_design="", aggregation_list=[]):
        if isinstance(data, list):
            seq_df = pd.DataFrame(data=data, columns=[SEQUENCES_COLUMN])
            self.original_sequences = data

        if isinstance(data, str):
            seq_df = pd.read_csv(data)


        if isinstance(data, pd.DataFrame):
            #print(data)
            seq_df = data

        self.library_design = library_design
        self.seq_df = seq_df

        # If the csv contains a count column there is no need re recount otherwise preform the count.
        if COUNT_COLUMN not in self.seq_df.columns:
            # Aggeregate similar sequences together
            aggregeate_by = [SEQUENCES_COLUMN] + aggregation_list
            self.seq_df = self.seq_df.groupby(aggregeate_by).size().to_frame(COUNT_COLUMN).reset_index()
        #print(self.seq_df)

    def preprocess_reads(self, preprocessor):
        """
        Preprocessor the library reads with a certain preprocess procedure defined by the given object. 
        :param preprocessor: The object in charge of the preproccessing. Should have a process method which takes a 
        string and return a string or None if the read should be removed from the library.
        :return: No return.
        """

        # Process sequences
        print("Number of reads before preprocessing : {}".format(self.seq_df[COUNT_COLUMN].sum()))
        self.seq_df[SEQUENCES_COLUMN] = \
            self.seq_df.apply(lambda row: preprocessor.process(row[SEQUENCES_COLUMN]), axis=1)

        # Remove nan values from the sequences column
        self.seq_df = self.seq_df.dropna(subset=[SEQUENCES_COLUMN])
        #self.seq_df = self.seq_df.sample(frac=0.20)
        print("Number of reads after preprcoessing : {}".format(self.seq_df[COUNT_COLUMN].sum()))

    # TODO I'm not sure this is the right place for this.
    # TODO Strong couplling between the reads and the matching might be better for an matching object to get this class.
    # TODO For now it stays!x
    def match_sequences(self, matcher):
        """"
        This method matches each sequence in the library with an oligo using the supplied matcher.
        The matcher should implement a find_best_oligo_match function.
        """
        print(self.seq_df)
        self.seq_df['variant_id'] = self.seq_df.swifter.apply(lambda row:  matcher.find_best_variant_match(row[SEQUENCES_COLUMN]), axis=1)

        # Sorting the data by the variant id.
        self.seq_df = self.seq_df.sort_values(ID_COLUMN).reset_index(drop=True)

    def get_unaggeregate_dataframe(self, df=None):
        '''
        Given a dataframe with 3 columns : SEQUENCE_COLUMN = ID_COLUMN = Count
        return a new dataframe where each row in the given dataframe is duplicated count times.
        :param df: The Dataframe you want to unaggergate
        :return: A dataframe where each row in the given dataframe is duplicated count times.
        '''
        if df is None:
            df = self.seq_df

        sequences = pd.np.empty(shape=df['count'].sum(), dtype=pd.np.object)
        ids = pd.np.empty(len(sequences), pd.np.int)

        agg_index = 0
        for index, row in df.iterrows():
            seq = row[SEQUENCES_COLUMN]
            c = row[COUNT_COLUMN]
            v_id = row[ID_COLUMN]

            sequences[agg_index:agg_index + c] = seq
            ids[agg_index:agg_index + c] = int(v_id)

            agg_index += c

        return_df = pd.DataFrame(pd.np.asarray([sequences, ids]).T, columns=[SEQUENCES_COLUMN, ID_COLUMN])
        return return_df

    def align_reads(self, library_design):
        '''
        Aligns the library reads to their matched variants.
        At the end of the run each read will have a cigar path to his matched variant.
        :param library_design:
        :return:
        '''
        if self.did_matching() is False:
            print("Can not align an un matched library")
            return

        def generate_cigar_path(row):
            variant_id = row[ID_COLUMN]
            if variant_id == -1:
                return -1

            read_sequence = row[SEQUENCES_COLUMN]
            variant_sequence = library_design.get_variant_sequence(variant_id).upper()

            # Aligning the read and the variant. variant=query, read=target
            index=read_sequence.find("AATTGAATGCTTGCTTGCCG")
            if index != -1:
                #align = edlib.align(read_sequence[20:index], variant_sequence[20:-20], task="path")
                align = edlib.align(read_sequence, variant_sequence, task="path")

            else:
                align = edlib.align(read_sequence, variant_sequence, task="path")
                #align = edlib.align(read_sequence[20:-20], variant_sequence[20:-20], task="path")

            # TODO find a better name.
            query_target_path = biology.parse_cigar(align['cigar'])
            return query_target_path

        self.seq_df['cigar_path'] = self.seq_df.apply(generate_cigar_path, axis=1)

    def save_library_state(self, name="library_reads.csv"):
        self.seq_df.to_csv("{}".format(name))

    # TODO : *HACK* Added an unaggregate flag to allow backward computability. Want to get rid of this ASAP.
    def get_unmatched_reads(self, unaggregate=False):
        unmatched_df = self.seq_df.loc[self.seq_df[ID_COLUMN] == -1]

        if unaggregate:
            unmatched_df = self.get_unaggeregate_dataframe(unmatched_df)

        return [Read(row) for _, row in unmatched_df.iterrows()]

    def get_unmatched_reads_count(self):
        """
        :return: The number of unmatched reads in the library. 
        """
        unmatched = self.get_unmatched_dataframe()
        return unmatched[COUNT_COLUMN].sum()

    def get_matched_reads_count(self):
        """
        :return: The number of matched reads in the library.
        """
        matched = self.get_matched_dataframe()
        return matched[COUNT_COLUMN].sum()

    # TODO : *HACK* Added an unaggregate flag to allow backward computability. Want to get rid of this ASAP.
    def get_matched_reads(self, unaggregate=False):
        """
        Get all the matched sequences in the library.
        :return: A generator of the sequences (iterator object)
        """
        matched_df = self.get_matched_dataframe()
        if unaggregate:
            matched_df = self.get_unaggeregate_dataframe(matched_df)

        return [Read(row) for _, row in matched_df.iterrows()]

    def iterate_matched_reads(self):
        matched_df = self.get_matched_dataframe()
        for _, row in matched_df.iterrows():
            yield Read(row)

    @staticmethod
    def get_weighted_count(df, pattern):
        base_count = df[SEQUENCES_COLUMN].str.count(pattern)
        base_count_weighted = base_count * df[COUNT_COLUMN]
        return base_count_weighted.sum()

    def get_total_matched_base_count(self):
        """
        Return the total number of bases among the matched reads.
        If we had 2 matched reads:"ACGT", "ACCGT" we return 9.
        :return: 
        """
        matched_df = self.get_matched_dataframe()

        return self.get_weighted_count(matched_df, "[A-Za-z]")

    def get_base_count_in_reads(self, base='A', matched=True):
        '''
        Get the number of occurrences of the specified base in the matched reads.
        Simply iterates over the reads and increment a counter every time the base was encountered
        :param base: (str) The base you want to count.
        :param matched : (bool) Get the count from the matched if True, get counts from all reads if False.
        :return: (int) The specified base count in the reads.
        '''
        pattern = "[{}]".format(base)
        if matched:
            matched_df = self.get_matched_dataframe()
            return self.get_weighted_count(matched_df, pattern)

        else:
            return self.get_weighted_count(self.seq_df, pattern)

    def get_base_count_in_design(self, base='A'):
        '''
        Get the number of occurrences of the specified base in the design. If the library
        was synthesized perfectly then get_base_count_in_design == get_base_count_in_reads for every base.
        Iterates over the variants and increment a counter by the number of reads
        for that variant every time the base was encountered.
        :param base: The base you want to count.
        :return: (int) The specified base count in .
        '''
        variants_df = self.library_design.get_df_copy()

        # Get base count by variants.
        base_count = variants_df['sequence'].str.count("[{}]".format(base))

        # Get variant counts in matched reads.
        matched_df = self.get_matched_dataframe()
        ids_count = matched_df.groupby(ID_COLUMN)[COUNT_COLUMN].sum()

        # Create count vector for multiplication with base count.
        count_vector = pd.np.zeros(len(base_count))
        count_vector[ids_count.index.values] = ids_count.values

        # Multiply base count by variant count.
        return (base_count * count_vector).sum()

    '''
        variants_df = self.library_design.get_df_copy()

        # Get base count by variants.
        base_count = variants_df['sequence'].str.count("[{}]".format(base))

        # Get variant counts in matched reads.
        matched_df = self.get_matched_dataframe()
        reads_grouped_df = matched_df.groupby(['variant_id']).size().to_frame('size')

        # Create count vector for multiplication with base count.
        count_vector = pd.np.zeros(len(base_count))
        count_vector[reads_grouped_df.index] = reads_grouped_df['size']

        # Multiply base count by variant count.
        return (base_count * count_vector).sum()
    '''

    def did_matching(self):
        if 'variant_id' in self.seq_df.columns:
            return True

        return False

    def did_edit_distance(self):
        if 'cigar_path' in self.seq_df.columns:
            return True

        return False

    def set_config_file(self, file_name):
        with open(file_name) as f:
            self.config = json.load(f)

    def get_config_file(self):
        if self.config:
            return self.config

        return None

    def set_original_sequences(self, sequences):
        print("This function is deprecated, original sequences are now resolved from existing data.")
        self.original_sequences = sequences.copy()

    def get_original_sequences(self):
        return self.seq_df[SEQUENCES_COLUMN].copy()

    def get_df_copy(self):
        return self.seq_df.copy()

    def __len__(self):
        return self.seq_df.shape(0)

    def __iter__(self):
        for _, read in self.seq_df.iterrows():
            yield Read(read)

    # ***************************************************************** Private methods *******************************
    def get_matched_dataframe(self):
        return self.seq_df.loc[self.seq_df[ID_COLUMN] != -1]

    def get_unmatched_dataframe(self):
        return self.seq_df.loc[self.seq_df[ID_COLUMN] == -1]
