""""
This class is intended to match a sequence with the most suitable variant.
"""
import time
import numpy as np

from src.library_design.library_design import LibraryDesign
from src.library_design import library_design_loader
from src.library_reads import library_reads_loader

NO_MATCH = -1


class VariantMatcher():
    def __init__(self, library_design):
        self.library_design = library_design
        # I wanted to us the oligo design class iterator and iterate over the class instead of
        # this explicit call but this is 25% faster.
        self.oligo_design_sequences = self.library_design.get_oligo_sequences()

        # If you want to add more matching algorithms simply write a function and add it to
        # this array using the API function add_aligner
        self.aligners = []

    def find_best_variant_match(self, read):
        '''
        Iterate over the matching functions in the matchers array. If a match was found
        it will be returned, other wise it will move on to the next matching function.
        If no match was found for all the function returns -1.
        :param sequence: The sequence we want to find the matching oligo for. 
        :return: the oligo index if one was found or -1 if no oligo match was found.
        '''
        indices = []
        for aligner in self.aligners:
            index = aligner.align_variant_to_read(indices, read)

            # Turning the index into
            if isinstance(index, int):
                index = [index]

            if len(index) > 1:
                indices = index
            elif index[0] != NO_MATCH:
                return index[0]

        return NO_MATCH

    def add_aligner(self, aligner):
        '''
        Add an aligner to
        :param aligner:
        :return:
        '''
        self.aligners.append(aligner)

    def full_match(self, sequence):
        '''
        Looks for a full match between the sequnces and the oligo design sequences.
        :param sequence: The sequence we want to find the matching oligo for.  
        :return: the oligo index if one was found or -1 if no oligo match was found.
        '''
        for index, oligo in enumerate(self.oligo_design_sequences):
            if oligo[35:] in sequence[35:]:
                return index
        return NO_MATCH


if __name__ == "__main__":
    print("\n\t##-- Running library_matcher.py --==\n")

    oligo_df = library_design_loader.load_oligo_design("../data/Twist_zohar_design.csv")
    oligo_design = LibraryDesign(oligo_df)
    oligo_matcher = VariantMatcher(oligo_design)
    sequences = library_loader.load_library("../data/merged", as_df=False)

    print(" - Data was loaded")
    t_start = time.time()

    for s in sequences:
        oligo_matcher.find_best_variant_match(s)

    print("Excecution time = {}".format(time.time() - t_start))

    # oligo_matcher.find_best_oligo_match("CCAATATCCTTAGCTGATCACCCATATGCCCACTAaaacaaacaaAGAGCGAATGGCGTAGTGCCGCATGAGGATTTCCTTATATTGGTGGTTAGTACGCATGCAATTAAAAATGTGCCCGGGCCACAGTTGACATTAGATATGG")