""""
This class represents the oligo design.
"""

import numbers

import src.library_design.library_design_loader
from src.library_design.variant import Variant

SEQUENCE_COLUMN = 'sequence'
BARCODE_COLUMN = 'barcode'
class LibraryDesign():

    # TODO might want to accept a sequqnce of oligos instead of a dataframe.
    def __init__(self, design_df):
        self.design_df = design_df

        # Name of columns will be in lower case for consistency.
        self.design_df.columns = map(str.lower, self.design_df.columns)

        # Uppercasing the barcode and sequence columns.
        self.design_df[SEQUENCE_COLUMN] = self.design_df[SEQUENCE_COLUMN].str.upper()
        self.design_df[BARCODE_COLUMN] = self.design_df[BARCODE_COLUMN].str.upper()

        # We will keep the sequences as a simple array for better iteration performance
        self.variants_seq = [x.upper() for x in design_df['sequence'].values]

    def get_variant_sequence(self, index, upper=True):
        variant_row = self.get_row_by_index(index)
        if upper:
            return variant_row['sequence'].upper()
        else:
            return variant_row['sequence']

    def get_all_variants_strings(self):
        variants = self.get_df_copy()[SEQUENCE_COLUMN]
        return variants.values

    def get_variant(self, id):
        variant_row = self.get_row_by_index(id)
        return Variant(variant_row)

    def get_oligo_sequences(self):
        return self.variants_seq

    def get_variant_barcode(self, variant):
        if isinstance(variant, numbers.Number):
            if variant > len(self):
                print ("-- Id {} is not a variant in the design".format(variant))
                return -1

            return self.design_df.iloc[variant - 1]['barcode']

        if isinstance(variant, basestring):
            index = -1
            for idx, seq in enumerate(self.variants_seq):
                if seq == variant.upper():
                    index = idx
                    break

            if index == -1:
                print("-- given sequence is not a variant in the design.")
                return -1
            else:
                row = self.design_df.iloc[index]
                return row['barcode']

    def has_barcode(self):
        if 'barcode' in self.design_df.columns:
            return True

        return False

    def get_library_barcodes(self, candidates=[]):
        if self.has_barcode():
            if not candidates:
                return self.design_df['barcode'].values
            else:
                return self.design_df.iloc[candidates]['barcode'].values
        else:
            print ("Requested barcode array for a library design with no barcode")
            return []

    def get_variant_by_id(self, id):
        variant_row = self.get_row_by_index(id)
        return Variant(variant_row)

    def get_variant_by_read(self, read):
        variant_id = read.get_variant_id()
        return self.get_variant_by_id(variant_id)

    def get_variants_strings_by_ids(self, ids):
        variants_str = self.design_df[SEQUENCE_COLUMN].copy()
        return variants_str[ids].values

    def get_barcode_strings_by_ids(self, ids):
        barcode_str = self.design_df[BARCODE_COLUMN].copy()
        return barcode_str[ids].values

    def get_number_of_variants(self):
        return len(self.design_df)

    def get_design_longest_sequence_len(self):
        max_length = self.design_df['sequence'].apply(lambda x: len(x)).max()
        return max_length

    def get_design_shortest_sequence_len(self):
        max_length = self.design_df['sequence'].apply(lambda x: len(x)).min()
        return max_length

    def get_design_sequences_average_len(self):
        max_length = self.design_df['sequence'].apply(lambda x: len(x)).mean()
        return max_length

    def get_df_copy(self):
        return self.design_df.copy()

    # +++++++++++++++++++++++++++++++++++++++++++++++++++ Private Methods ++++++++++++++++++++++++++++++++++++++++++++++
    def get_row_by_index(self, index):
        return self.design_df.iloc[index]

    def __len__(self):
        return len(self.variants_seq)

    # def __getitem__(self, index):
    #     return Variant(self.oligo_df.iloc[index])

    def __iter__(self):
        for index, variant_row in self.design_df.iterrows():
            yield Variant(variant_row)


def check_enumeration(design):
    for index, o in enumerate(design):
        print (o)

def check_barcode(design):
    if design.has_barcode():
        print("Variant design has Barcode!")
        print (design.get_variant_barcode(3))
        print (design.get_variant_barcode("CCAATATCCTTAGCTGATCACCCATATGCCCACTAaaacaaatggAGAGCGAATGGCGTAGTGCCGCAACGTTGTTCTCATCGTCGATAAAATGGCATGAGAGTTGCTGTGTTTTAGCTGCCCGGGCCACAGTTGACATTAGATATGG"))
    else:
        print("Variant design has no Barcode!")

def test_library_barcodes():
    ld = get_library_design_sample()
    try:
        barcodes = ld.get_library_barcodes()
    except:
        print ("Get library barcodes fail")
        raise

    try:
        barcodes_candidate = ld.get_library_barcodes([0])
    except:
        print ("Get library barcodes by candidates failed")
        raise

    if barcodes[0] != barcodes_candidate[0]:
        print ("password")


def get_library_design_sample(file="../data/Twist_zohar_design.csv"):
    design_df = src.library_design.library_design_loader.load_oligo_design(file)
    oligo_design = LibraryDesign(design_df)

    return oligo_design

if __name__ == "__main__":
    print("\n\t == Testing Oligo design class ==\n")
    design_df = src.library_design.library_design_loader.load_oligo_design("../data/Twist_zohar_design.csv")

    oligo_design = LibraryDesign(design_df)
    check_barcode(oligo_design)

    print("\n\t == Testing Oligo design class Done ==")
