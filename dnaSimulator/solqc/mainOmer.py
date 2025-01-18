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
from src.library_reads import library_loader
from src.library_reads.library_reads import LibraryReads
from src.library_reads.preprocessor import FixPreProcessor
from src.utils.pdf_generator import PDFGenerator
import src.config as config

# Plotting setup
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

sns.set_palette('colorblind')
plt.rc("axes.spines", top=False, right=False) # Remove top and right borders. Looks nicer.

mpl.rcParams['savefig.dpi'] = 300 # Save plot will be in a higher resolution.


# TODO move file names to another file.


#DATASET_FOLDER = "/Users/omersabary/Desktop/LibAnalyzer/Library-Analyzer-master/data/Erlich_unfiltered/"
#DATASET_FOLDER = "/Users/omersabary/Desktop/LibAnalyzer/Library-Analyzer-master/data/Erlich/"
#DATASET_FOLDER = "/Users/omersabary/Desktop/LibAnalyzer/Library-Analyzer-master/data/Erlich_10/"

#DATASET_FOLDER = "/Users/omersabary/Desktop/Luis/"
#DATASET_FOLDER= "/Users/omersabary/Documents/Luis_Matching/"
DATASET_FOLDER = "/Users/omersabary/Dropbox/Grass_data/TrellisBMAdata/"
DATASET_FOLDER = "/Users/omersabary/Dropbox/SOLQC_Analysis_new_perfect/Grass/"
DATASET_FOLDER = "/Users/omersabary/Dropbox/Grass_data/split_data_sets/Pfitser/"
DATASET_FOLDER = "/Users/omersabary/Downloads/solqc_verifica/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/dataPilotPool/DataForSOLQC_frac0.02/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/dataPilotPool/simulatorValdiation/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/nanoporedataMay10/solqc_files/"
#DATASET_FOLDER = "/Users/omersabary/Dropbox/Grass_data/SOLQC_Grass_Omer/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/nanoporedataMay10/solqc_files/SOLQC_file_primers_10_full_rev_com/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/nanoporedataMay10/solqc_files/SOLQC_file_primers_10_full_rev_com/psuedo_10_rev_com/"
DATASET_FOLDER = "/Users/omersabary/Desktop/DaniellaDvir/simuComparison_Omer/original_test_1/"

#DATASET_FOLDER = ""
#DATASET_FOLDER = "/Users/omersabary/Documents/HosseinData/Files/hossein_final/"
#DATASET_FOLDER = "/Users/omersabary/Documents/HosseinData/Files/"

#DATASET_FOLDER = "/Users/omersabary/Documents/Reinhart/"
#DATASET_FOLDER = "/Users/omersabary/Documents/Reinhart_unfiltered/"


#DATASET_FILE = ""
#DATASET_FILE = "erlich"
#DATASET_FILE = "hossein"
#DATASET_FILE = "reinhart"
DATASET_FILE = "pfitser"
DATASET_FILE = "evyatSOLQC"
DATASET_FILE = "dvir"
#DATASET_FILE = "reinhard"



DATASET_READS_SUFFIX = "_reads.fastq"
DATASET_DESIGN_SUFFIX = "_design_n.csv"
DATASET_DESIGN_SUFFIX = "_design.csv"
DATASET_CONFIG_SUFFIX = "_config.json"


# A const to decide whether to use existing files or re calculate. [Matching and Alignment].
OVERRIDE = False
EDIT_DISTANCE = True
# An int to allow a smaller chunk of reads from the fastq file.
NREADS = -1


def get_library_reads_and_design(reads_file, design_file, override=False):
    print("\n| Setting up library design and reads object:")
    library_design_df = library_design_loader.load_oligo_design(design_file)
    library_design = LibraryDesign(library_design_df)
    print(" - Generated the design obj")

    read_file_noprefix = reads_file[:-12]
    print(read_file_noprefix)
    print(reads_file)
    aligned_file = Path("{}_after_alignment.csv".format(read_file_noprefix))
    #aligned_file = Path("{}_reads_matching.csv".format(read_file_noprefix))
    matched_file = Path("{}_after_matching.csv".format(read_file_noprefix))
    matched_file = Path("{}_after_matching.csv".format(read_file_noprefix))
    config_file = Path(read_file_noprefix + DATASET_CONFIG_SUFFIX)
    print(config_file)

    sequences = library_loader.load_library(DATASET_FOLDER+ DATASET_FILE+DATASET_READS_SUFFIX, True, NREADS)
    #print(sequences)
    #exit(0)
    #sequences = sequences.sample(frac=0.02)
    #print(sequences)
    #exit(0)
    if aligned_file.exists() and not override:
        print("aligned")
        library_reads = LibraryReads(str(aligned_file), library_design, aggregation_list=['variant_id'])
    elif matched_file.exists() and not override:
        print("matched")
        library_reads = LibraryReads(str(matched_file), library_design, aggregation_list=['variant_id'])
    else:
        library_reads = LibraryReads(sequences, library_design)

    if config_file.exists():
        print("start")
        config.set_configuration_file(config_file)
        print(" - Config file set")
    else:
        assert 'Config file does not exist'

    print(" - Generated the reads obj")
    return library_reads, library_design


def get_library_reads_and_design_simple(design_file, read_metafile):
    print(read_metafile)
    library_design_df = library_design_loader.load_oligo_design(design_file)
    library_design = LibraryDesign(library_design_df)
    print(" - Library design initialized successfully.")

    reads_file_names = get_filenames_from_file(read_metafile)
    sequences = library_loader.load_library(reads_file_names)

    library_reads = LibraryReads(sequences, library_design)
    print(" - Library reads initialized successfully.")

    return library_reads, library_design


def get_filenames_from_file(file_name):
    with open(file_name, 'r') as file:
        files_names = file.read().splitlines()

    return files_names


def preprocess_reads(library_reads):
    print(" - Pre Processing")
    processor = FixPreProcessor(config.get_configuration())
    library_reads.preprocess_reads(processor)


def match_reads_to_design(library_reads, library_design, aligners_names=None):
    if aligners_names is None:
        aligners_names = ["BarcodeAligner"]
        aligners_names = ["EditDistanceAligner"]
    
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


def align_reads_to_variants(library_reads, library_design):
    print(" - Aligning")
    library_reads.align_reads(library_design)

    library_reads.save_library_state("{}_after_alignment.csv".format(DATASET_FOLDER+DATASET_FILE))
    print(" - - Finished read alignment")


def analyze_library(library_reads, library_design, analyzers_array=[], report_id=""):
    contents = []
    analyzers = AnalyzerFactory.create_analyzers(analyzers_array)
    #analyzers.sort() # sort analyzers by their display_rank.

    for analyzer in analyzers:
        contents.append(analyzer.analyze(library_reads, library_design))
        print(" - - Finished with : {}".format(analyzer.name))

    print(analyzers_array)

    save_path = config.get_deliverable_path()
    with PDFGenerator("{}/report.pdf".format(save_path)) as pdf_gen:
        for content in contents:
            if content:
                pdf_gen.add_content_array(content)


# TODO : this is for working with CONST values, but we need to eliminate this.
def start_analyzer_old():
    lr_obj, ld_obj = get_library_reads_and_design(DATASET_FOLDER + DATASET_FILE + DATASET_READS_SUFFIX,
                                                  DATASET_FOLDER + DATASET_FILE + DATASET_DESIGN_SUFFIX,
                                                  OVERRIDE)
    config.set_configuration_file(DATASET_FOLDER + DATASET_FILE+ DATASET_CONFIG_SUFFIX)
    if not lr_obj.did_matching() or OVERRIDE:
        print("matching was not done")
        preprocess_reads(lr_obj)
        match_reads_to_design(lr_obj, ld_obj)

    #if EDIT_DISTANCE and (not lr_obj.did_edit_distance() or OVERRIDE):
    print("omer")
    align_reads_to_variants(lr_obj, ld_obj)

    # This delay fix a strange error. Probably something with reading from a file we just finished writing to.
    time.sleep(2)
    analyzers_array=[AnalyzersNames.GENERAL_ANALYZER,
        AnalyzersNames.BASE_DISTRIBUTION_PER_POSITION,
        AnalyzersNames.VARIANT_DISTRIBUTION_ANALYZER,
        #AnalyzersNames.VARIANT_DISTRIBUTION_BY_GC_ANALYZER,
        AnalyzersNames.READS_LENGTH_ANALYZER,
        AnalyzersNames.ERROR_RATES_ANALYZER,
        AnalyzersNames.BASE_DEPENDENT_ANALYZER,
        AnalyzersNames.ERROR_PER_POSITION_ANALYZER]
        #AnalyzersNames.LONG_DELETION_ANALYZER,
        #AnalyzersNames.CUMULATIVE_ERROR_ANALYZER,
        #AnalyzersNames.ERROR_RATE_BY_GC_ANALYZER]
        #AnalyzersNames.LONG_DELETION_ANALYZER]
    analyze_library(lr_obj, ld_obj, analyzers_array=analyzers_array)


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
    lr_obj, ld_obj = get_library_reads_and_design_simple(args.design, args.reads)
    # Preforming variants <-> reads matching.
    if not lr_obj.did_matching() or OVERRIDE:
        preprocess_reads(lr_obj)
        match_reads_to_design(lr_obj, ld_obj, args.aligners)

    # TODO : We should have a mechansim which enables the use of previous run (similar to the one we used in debugging)
    # Preforming variants <-> reads alignment.
    if should_align(args.analyzers):
        align_reads_to_variants(lr_obj, ld_obj)

    # This delay fix a strange error.
    time.sleep(2)

    # Analyzing the gathered data of the matching and alignment.
    analyze_library(lr_obj, ld_obj, analyzers_array=args.analyzers, report_id=args.id)


if __name__ == "__main__":
    print("Starting QC Analysis")
    t_start = time.time()
    #start_analyzer(sys.argv[1:])
    start_analyzer_old()
    print("\n -(\) - Finished QC analysis in : {:.2f} seconds.".format(time.time() - t_start))
    # end of main

# To analyze the library we need both the design and the library reads object.
# The design can be loaded regularly.
# What we need is the variant_id and the cigar path.
