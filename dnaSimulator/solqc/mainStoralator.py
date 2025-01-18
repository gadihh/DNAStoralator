import time
from pathlib import Path
import sys
import os
import datetime

from src.aligners.aligner_factory import AlignerFactory
from src.analyzers.analyzer import AnalyzersNames
from src.analyzers.analyzer_factory import AnalyzerFactory
from src.utils.qc_argumnet_parser import parse_arguments
from src.aligners.library_matcher import VariantMatcher
from src.library_design import library_design_loader
from src.library_design.library_design import LibraryDesign
from src.library_reads import library_reads_loader
from src.library_reads.library_reads import LibraryReads
from src.library_reads.preprocessor import FixPreProcessor
from src.utils.pdf_generator import PDFGenerator
import src.config as config

# Plotting setup
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

sns.set_palette('muted')
plt.rc("axes.spines", top=False, right=False) # Remove top and right borders. Looks nicer.
mpl.rcParams['savefig.dpi'] = 300 # Save plot will be in a higher resolution.


def preprocess_reads(library_reads):
    print(" - Pre Processing")
    processor = FixPreProcessor(config.get_configuration())
    library_reads.preprocess_reads(processor)


def load_design(design_input):
    library_design_df = library_design_loader.load_oligo_design(design_input)
    library_design = LibraryDesign(library_design_df)
    print(" - Library design initialized successfully.")

    return library_design


def load_reads(reads_input, library_design):
    reads_init_input = library_reads_loader.load_library(reads_input)
    library_reads = LibraryReads(reads_init_input, library_design)
    print(" - Library reads initialized successfully.")

    return library_reads


def match_reads_to_design(library_reads, library_design, aligners_names=None):
    if aligners_names is None:
        aligners_names = ["BarcodeAligner"]

    print(" - Matching")
    variant_matcher = VariantMatcher(library_design)

    aligners = AlignerFactory.create_aligners(aligners_names, library_design)
    for aligner in aligners:
        variant_matcher.add_aligner(aligner)

    library_reads.match_sequences(variant_matcher)
    # library_reads.save_library_state("{}_after_matching.csv".format(DATASET_FOLDER+DATASET_FILE))

    matched_count = library_reads.get_matched_reads_count()
    un_matched_count = library_reads.get_unmatched_reads_count()
    print(" - - managed to match - {}".format(matched_count))
    print(" - - % of successful match - {}".format(matched_count / (matched_count + un_matched_count)))


def should_align(analyzers_names):
    '''
    Given an array of analyzers checks whether any of the analyzers requires alignment.
    :param analyzers_names: The analyzers to check
    :return: True if any of the analyzers requires alignment false otherwise.
    '''
    analyzers = AnalyzerFactory.create_analyzers(analyzers_names)
    for analyzer in analyzers:
        if analyzer.requires_alignment:
            return True

    return False


def align_reads_to_variants(library_reads, library_design):
    print(" - Aligning")
    library_reads.align_reads(library_design)
    print(" - - Finished read alignment")


def analyze_library(library_reads, library_design, analyzers_array=[], report_id=""):
    contents = []
    analyzers = AnalyzerFactory.create_analyzers(analyzers_array)
    analyzers.sort() # sort analyzers by their display_rank.

    for analyzer in analyzers:
        contents.append(analyzer.analyze(library_reads, library_design))
        print(" - - Finished with : {}".format(analyzer.name))

    print(analyzers_array)

    save_path = config.get_deliverable_path()
    with PDFGenerator("{}/report.pdf".format(save_path)) as pdf_gen:
        for content in contents:
            if content:
                pdf_gen.add_content_array(content)

def analyze_library_storalator(library_reads, library_design, analyzers_array=[], report_id=""):
    contents = []
    analyzers = AnalyzerFactory.create_analyzers(analyzers_array)
    #analyzers.sort() # sort analyzers by their display_rank.


    for analyzer in analyzers:
        #dickeys, dicvals = (analyzer.analyze(library_reads, library_design))
        #print(dickeys)
        #print(dicvals)

        contents.append(analyzer.analyze(library_reads, library_design))
        print(" - - Finished with : {}".format(analyzer.name))
    #print(library_design)
    #contents.append(analyzers[2].analyze(library_reads, library_design))
    #print(" - - Finished with : {}".format(analyzer.name))

    print(analyzers_array)

    save_path = config.get_deliverable_path()


    #with PDFGenerator("{}/report.pdf".format(save_path)) as pdf_gen:
    #    for content in contents:
    #        if content:
    #            pdf_gen.add_content_array(content)


def analyze_matched_library():
    pass


def start_analyzer(arguments):
    args = parse_arguments(arguments)

    # Setting a global configuration file.
    config.set_configuration_file(args.config)
    config.set_deliverable_path(args.time)

    # Temporary images will be saved in the temp folder, if not temp folder exists we create one.
    if not os.path.isdir("temp"):
        os.makedirs("temp")

    # Deliverables Will be placed into the deliverable folder.
    if not os.path.isdir("deliverable"):
        os.makedirs("deliverable")

    # Create the deliverable folder for the current run.
    os.makedirs("{}".format(config.get_deliverable_path()))

    # Building the library reads and library design objects.
    library_design = load_design(design_input=args.design)
    library_reads = load_reads(reads_input=args.reads, library_design=library_design)
    # Preforming variants <-> reads matching.
    if not library_reads.did_matching():
        preprocess_reads(library_reads)
        match_reads_to_design(library_reads, library_design, args.aligners)

    # Preforming variants <-> reads alignment.
    if should_align(args.analyzers) and not library_reads.did_edit_distance():
        align_reads_to_variants(library_reads, library_design)

    # This delay fix a strange error.
    time.sleep(2)

    # Analyzing the gathered data of the matching and alignment.
    analyze_library(library_reads, library_design, analyzers_array=args.analyzers, report_id=args.id)

def start_analyzer_storalator(arguments):
    print("omer")
    print(arguments)
    args = parse_arguments(arguments)
    print("omer")

    # Setting a global configuration file.
    config.set_configuration_file(args.config)
    config.set_deliverable_path(args.time)

    # Temporary images will be saved in the temp folder, if not temp folder exists we create one.
    if not os.path.isdir("temp"):
        os.makedirs("temp")

    # Deliverables Will be placed into the deliverable folder.
    if not os.path.isdir("deliverable"):
        os.makedirs("deliverable")

    # Create the deliverable folder for the current run.
    os.makedirs("{}".format(config.get_deliverable_path()))

    # Building the library reads and library design objects.
    library_design = load_design(design_input=args.design)
    print(args.reads)
    library_reads = load_reads(reads_input=args.reads, library_design=library_design)
    # Preforming variants <-> reads matching.
    if not library_reads.did_matching():
        preprocess_reads(library_reads)
        match_reads_to_design(library_reads, library_design, args.aligners)

    # Preforming variants <-> reads alignment.
    if should_align(args.analyzers) and not library_reads.did_edit_distance():
        align_reads_to_variants(library_reads, library_design)

    # This delay fix a strange error.
    time.sleep(2)

    # Analyzing the gathered data of the matching and alignment.
    analyze_library_storalator(library_reads, library_design, analyzers_array=[AnalyzersNames.ERROR_RATES_ANALYZER_STORALATOR, AnalyzersNames.BASE_DEPENDENT_ANALYZER_STORALATOR, AnalyzersNames.VARIANT_DISTRIBUTION_ANALYZER_STORALATOR], report_id=args.id)



if __name__ == "__main__":
    print("Starting QC Analysis")
    t_start = time.time()
    start_analyzer_storalator(sys.argv[1:])
    print("\n -(\) - Finished QC analysis in : {:.2f} seconds.".format(time.time() - t_start))
    # end of main

# To analyze the library we need both the design and the library reads object.
# The design can be loaded regularly.
# What we need is the variant_id and the cigar path.
