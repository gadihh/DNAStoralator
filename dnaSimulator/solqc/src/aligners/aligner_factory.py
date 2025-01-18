import src.config as config
from src.aligners.aligner import AlignersNames
from src.aligners.barcode_aligner import BarcodeAligner
from src.aligners.edit_distance_aligner import EditDistanceAligner
from src.aligners.n_edit_barcode_aligner import NEditBarcodeAligner

class AlignerFactory(object):

    @staticmethod
    def create_aligners(aligners_names, library_design):
        configuration = config.get_configuration()
        aligners = []
        for name in aligners_names:
            if name == AlignersNames.BARCODE_ALIGNER:
                # TODO wanted keys should be in the aligner itself.
                wanted_keys = ['barcode_start', 'barcode_end']
                aligner = BarcodeAligner(library_design,
                                         dict((k, configuration[k]) for k in wanted_keys if k in configuration))
                aligners.append(aligner)

            if name == AlignersNames.N_EDIT_BARCODE_ALIGNER:
                wanted_keys = ['barcode_start', 'barcode_end', 'barcode_tolerance']
                aligner = NEditBarcodeAligner(library_design,
                                              dict((k, configuration[k]) for k in wanted_keys if k in configuration))
                aligners.append(aligner)

            if name == AlignersNames.EDIT_DISTANCE_ALIGNER:
                aligner = EditDistanceAligner(library_design, {})
                aligners.append(aligner)

        return aligners
