from progress.bar import Bar


class AnalyzersNames:

    NEW_BASE_DEPENDENT_ANALYZER         = "NewBaseDependentAnalyzer"
    NEW_CUMMULATIVE                     = "NewCumulativeErrorAnalyzer"
    BASE_DEPENDENT_ANALYZER             = "BaseDependentAnalyzer"
    BASE_DISTRIBUTION_PER_POSITION      = "BaseDistributionPerPositionAnalyzer"
    CUMULATIVE_ERROR_ANALYZER           = "CumulativeErrorAnalyzer"
    DELETION_ANALYZER                   = "DeletionAnalyzer"
    DELETION_POSITION_PER_BASE_ANALYZER = "DeletionPositionPerBaseAnalyzer"

    ERROR_PER_POSITION_ANALYZER         = "ErrorPerPositionAnalyzer"
    ERROR_RATE_BY_GC_ANALYZER           = "ErrorRateByGCAnalyzer"
    ERROR_RATE_PER_BASE_TABLE_ANALYZER  = "ErrorRatePerBaseTableAnalyzer"
    ERROR_RATES_ANALYZER                = "ErrorRatesAnalyzer"
    FREQUENCY_ANALYZER                  = "FrequencyAnalyzer"

    GENERAL_ANALYZER                    = "GeneralAnalyzer"
    LONG_DELETION_ANALYZER              = "LongDeletionAnalyzer"
    MATCHING_ANALYZER                   = "MatchingAnalyzer"
    READS_LENGTH_ANALYZER               = "ReadsLengthAnalyzer"
    VARIANT_DISTRIBUTION_BY_GC_ANALYZER = "VariantDistributionByGCAnalyzer"

    VARIANT_DISTRIBUTION_ANALYZER       = "VariantDistributionAnalyzer"
    
    ERROR_RATES_ANALYZER_STORALATOR     = "ErrorRatesAnalyzer_storalator"

    BASE_DEPENDENT_ANALYZER_STORALATOR = "NewBaseDependentAnalyzer_Storalator"

    VARIANT_DISTRIBUTION_ANALYZER_STORALATOR = "VariantDistributionAnalyzer_Storalator"


class Analyzer(object):
    name = "Analyzer"
    display_rank = 1000
    requires_alignment = True

    def __init__(self):
        pass

    def __lt__(self, other):
        return self.display_rank < other.display_rank

    def analyze(self, library_reads, library_design):
        """
        Analyzes a certain criteria of the library
        :param library_reads: A Library reads object.
        :return: A content object which can be used to visualize the analysis.
        """
        print("This method should be implemented by the analyzer object.")
        return False

    def get_progress_bar(self, length):
        return Bar('{} Analyzing'.format(str(self)), max=length)

    def __str__(self):
        return "Analyzer"
