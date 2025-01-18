import numpy as np

from src.aligners.aligner import Aligner


class BarcodeAligner(Aligner):
    def __init__(self, library_design, data):
        """
        The Barcode Aligner class uses an exact barcode equality to preform a matching between
        A Read and a Varaint.
        :param library_design: A LibraryDesign object.
        :param data: A dictionary with the relevant data for the aligner : [barcode_start, barcode_end].
        """
        super(BarcodeAligner, self).__init__()
        if not library_design.has_barcode():
            print ("The library design sent to barcode aligner has no barcode.")

        self.library_design = library_design
        self.data = data

        barcode_start = data['barcode_start']
        barcode_end = data['barcode_end']
        self.extract_barcode = lambda x: x[barcode_start:barcode_end]

        barcodes = library_design.get_df_copy()['barcode']
        self.barcodes_hash = {barcode:i for i, barcode in enumerate(barcodes)}

    def align_variant_to_read(self, variants_id, read):
        # If we have specific variants id we pick them of our barcodes hash, other wise use the original hash.
        if variants_id is not None and len(variants_id) > 0:
            barcodes_hash = {barcode:i for barcode, i in self.barcodes_hash.items() if i in variants_id}
        else:
            barcodes_hash = self.barcodes_hash

        # barcodes = library_design.get_barcode_strings_by_ids(variants_id)
        # barcodes_hash = {barcode:i for i, barcode in enumerate(barcodes)}

        reads = self.extract_barcode(read)
        if isinstance(reads, str):
            reads = [reads]

        for read in reads:
            val = barcodes_hash.get(read)
            if val is not None:
                return val

        return -1

    def get_library_design_barcodes(self, candidates=[]):
        return self.library_design.get_library_barcodes(candidates)


if __name__ == "__main__":
    pass
