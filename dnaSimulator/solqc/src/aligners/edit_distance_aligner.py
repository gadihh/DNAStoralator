'''

'''
# - Remove this remark once you've added the analyzer to the AnalyzerFactory.
from src.aligners.aligner import Aligner
import edlib
import src.config as config

class EditDistanceAligner(Aligner):
    name="Aligner"

    def __init__(self, library_design, data):
        """
        The Barcode Aligner class uses an exact barcode equality to preform a matching between
        A Read and a Varaint.
        :param library_design: A LibraryDesign object.
        :param data: A dictionary with the relevant data for the aligner : [barcode_start, barcode_end].
        """
        super(EditDistanceAligner, self).__init__()

        self.library_design = library_design
        self.data = data

        self.variants = library_design.get_all_variants_strings()

    def align_variant_to_read(self, variants_id, read):
        if variants_id is not None and len(variants_id) > 1:
            variants = self.variants[variants_id]
        else:
            variants = self.variants

        first_variant = variants[0]
        distance = edlib.align(read, first_variant)['editDistance']
        min_edit_d = distance
        min_index = 0

        for index, variant in enumerate(variants):
            distance = edlib.align(read, variant, k=min_edit_d)['editDistance']
            if -1 < distance < min_edit_d:
                min_edit_d = distance
                min_index = index

        return min_index