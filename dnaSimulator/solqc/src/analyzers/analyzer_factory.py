from src.analyzers.analyzer import AnalyzersNames
from src.analyzers.error_statistics.base_dependent_analyzer import BaseDependentAnalyzer
from src.analyzers.composition_statistics.base_distribution_per_position_analyzer import BaseDistributionPerPositionAnalyzer
from src.analyzers.error_statistics.cumulative_error_analyzer import CumulativeErrorAnalyzer
from src.analyzers.error_statistics.deletion_analyzer import DeletionAnalyzer
from src.analyzers.error_statistics.deletion_position_per_base_analyzer import DeletionPositionPerBaseAnalyzer
from src.analyzers.error_statistics.error_per_position_analyzer import ErrorPerPositionAnalyzer
from src.analyzers.error_statistics.error_rate_by_gc_analyzer import ErrorRateByGCAnalyzer
from src.analyzers.error_statistics.error_rate_per_base_table_analyzer import ErrorRatePerBaseTableAnalyzer
from src.analyzers.error_statistics.error_rates_analyzer import ErrorRatesAnalyzer
from src.analyzers.raw_data.frequency_analyzer import FrequencyAnalyzer
from src.analyzers.composition_statistics.general_analyzer import GeneralAnalyzer
from src.analyzers.error_statistics.long_deletion_analyzer import LongDeletionAnalyzer
from src.analyzers.raw_data.matching_analyzer import MatchingAnalyzer
from src.analyzers.composition_statistics.reads_length_analyzer import ReadsLengthAnalyzer
from src.analyzers.composition_statistics.variant_distribution_by_GC_analyzer import VariantDistributionByGCAnalyzer
from src.analyzers.composition_statistics.variant_distriution_analyzer import VariantDistributionAnalyzer
from src.analyzers.error_statistics.base_dependent_analyzer_new import NewBaseDependentAnalyzer
from src.analyzers.error_statistics.cummulative_error_analyzer_new import NewCumulativeErrorAnalyzer
from src.analyzers.error_statistics.error_rates_analyzer_storalator import ErrorRatesAnalyzer_storalator
from src.analyzers.error_statistics.base_dependent_analyzer_new_storalator import NewBaseDependentAnalyzer_Storalator
from src.analyzers.composition_statistics.variant_distriution_analyzer_storalator import VariantDistributionAnalyzer_Storalator

class AnalyzerFactory(object):
    @staticmethod
    def create_analyzers(analyzers_names):
        analyzers = []
        for name in analyzers_names:
            if name == AnalyzersNames.BASE_DEPENDENT_ANALYZER:
                analyzers.append(BaseDependentAnalyzer())

            if name == AnalyzersNames.BASE_DISTRIBUTION_PER_POSITION:
                analyzers.append(BaseDistributionPerPositionAnalyzer())

            if name == AnalyzersNames.CUMULATIVE_ERROR_ANALYZER:
                analyzers.append(CumulativeErrorAnalyzer())

            if name == AnalyzersNames.DELETION_ANALYZER:
                analyzers.append(DeletionAnalyzer())

            if name == AnalyzersNames.DELETION_POSITION_PER_BASE_ANALYZER:
                analyzers.append(DeletionPositionPerBaseAnalyzer())

            if name == AnalyzersNames.ERROR_PER_POSITION_ANALYZER:
                analyzers.append(ErrorPerPositionAnalyzer())

            if name == AnalyzersNames.ERROR_RATE_BY_GC_ANALYZER:
                analyzers.append(ErrorRateByGCAnalyzer())

            if name == AnalyzersNames.ERROR_RATE_PER_BASE_TABLE_ANALYZER:
                analyzers.append(ErrorRatePerBaseTableAnalyzer())

            if name == AnalyzersNames.ERROR_RATES_ANALYZER:
                analyzers.append(ErrorRatesAnalyzer())

            if name == AnalyzersNames.FREQUENCY_ANALYZER:
                analyzers.append(FrequencyAnalyzer())

            if name == AnalyzersNames.GENERAL_ANALYZER:
                analyzers.append(GeneralAnalyzer())

            if name == AnalyzersNames.LONG_DELETION_ANALYZER:
                analyzers.append(LongDeletionAnalyzer())

            if name == AnalyzersNames.MATCHING_ANALYZER:
                analyzers.append(MatchingAnalyzer())

            if name == AnalyzersNames.READS_LENGTH_ANALYZER:
                analyzers.append(ReadsLengthAnalyzer())

            if name == AnalyzersNames.VARIANT_DISTRIBUTION_BY_GC_ANALYZER:
                analyzers.append(VariantDistributionByGCAnalyzer())

            if name == AnalyzersNames.VARIANT_DISTRIBUTION_ANALYZER:
                analyzers.append(VariantDistributionAnalyzer())

            if name == AnalyzersNames.NEW_BASE_DEPENDENT_ANALYZER:
                analyzers.append(NewBaseDependentAnalyzer())

            if name == AnalyzersNames.NEW_CUMMULATIVE:
                analyzers.append(NewCumulativeErrorAnalyzer())

            if name == AnalyzersNames.ERROR_RATES_ANALYZER_STORALATOR:
                analyzers.append(ErrorRatesAnalyzer_storalator())
            
            if name == AnalyzersNames.BASE_DEPENDENT_ANALYZER_STORALATOR: 
                analyzers.append(NewBaseDependentAnalyzer_Storalator())

            if name == AnalyzersNames.VARIANT_DISTRIBUTION_ANALYZER_STORALATOR:
                analyzers.append(VariantDistributionAnalyzer_Storalator())
        return analyzers

if __name__ == '__main__':
    analyzers_list = [
        "BaseDependentAnalyzer",
        "BaseDistributionPerPositionAnalyzer",
        "CumulativeErrorAnalyzer",
        "DeletionAnalyzer",
        "DeletionPositionPerBaseAnalyzer",
        "ErrorPerPositionAnalyzer",
        "ErrorRateByGCAnalyzer",
        "ErrorRatePerBaseTableAnalyzer",
        "ErrorRatesAnalyzer",
        "FrequencyAnalyzer",
        "GeneralAnalyzer",
        "LongDeletionAnalyzer",
        "MatchingAnalyzer",
        "ReadsLengthAnalyzer",
        "VariantDistributionByGCAnalyzer",
        "VariantDistributionAnalyzer"
    ]

    af = AnalyzerFactory()
    af.create_analyzers(analyzers_list)