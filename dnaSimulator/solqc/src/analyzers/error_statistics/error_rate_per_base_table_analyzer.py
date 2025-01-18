"""
Error Table Analyzer
This class is in charge of analyzing the amount of errors (missmatch, deletion
insertion, long deletion and long insertion) that was found in the library and
display them in a tabular fashion.
"""

# Global Imports
import matplotlib.pyplot as plt

# Local Imports
from src.analyzers.analyzer import Analyzer
from src.utils.biology import DELETION, INSERTION, MISMATCH
from src.utils.content import Content


class DNAError ():
    n_missmatch = 0
    n_deletion = 0
    n_insertion = 0
    n_long_deletion = 0
    n_long_insertion = 0

    def increment_missmatch(self, by=1):
        self.n_missmatch += by

    def increment_n_deletion(self, by=1):
        self.n_deletion += by

    def increment_n_insertion(self, by=1):
        self.n_insertion += by

    def increment_n_long_deletion(self, by=1):
        self.n_long_deletion += by

    def increment_n_long_insertion(self, by=1):
        self.n_long_insertion += by

    def get_data(self):
        return (self.n_missmatch,
                self.n_deletion,
                self.n_insertion,
                self.n_long_deletion,
                self.n_long_insertion)


class ErrorRatePerBaseTableAnalyzer(Analyzer):
    name = "Error Table Analyzer"
    requires_alignment = True

    def __init__(self):
        super().__init__()

    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance():
            print("You can not perform deletion analyzing on reads that did not go through edit distance")
            return None

        self.error_table = {
            "A" : DNAError(),
            "C" : DNAError(),
            "G" : DNAError(),
            "T" : DNAError()
        }

        matched_reads = library_reads.get_matched_reads()

        for read in matched_reads:
            cigar = read.get_cigar_path()
            variant = library_design.get_variant_by_read(read)

            var_cigar = cigar.replace(INSERTION, "")
            read_cigar = cigar.replace(DELETION, "")
            read_count = read.get_row_count()


            self.update_read_errors(read(), read_cigar,read_count)
            self.update_var_errors(variant(), var_cigar, read_count)

        return self.generate_content(library_reads)

    def update_read_errors(self, read, read_cigar, count):
        long_insertion = False
        for i, s in enumerate(read_cigar):
            if s == INSERTION:
                if not long_insertion:
                    if i == len(read_cigar) - 1 or read_cigar[i + 1] != INSERTION:
                        self.error_table[read[i]].increment_n_insertion(count)
                    else:
                        self.error_table[read[i]].increment_n_long_insertion(count)
                else:  # if the long_deletion flag is on there is no updating.
                    pass
                # Set long flag deletion to on.
                    long_insertion = True
            else:
                long_insertion = False

    def update_var_errors(self, variant, variant_cigar, count):
        long_deletion = False
        for i, s in enumerate(variant_cigar):
            if s == DELETION:
                if not long_deletion:
                    if i == len(variant_cigar) - 1 or variant_cigar[i+1] != DELETION:
                        self.error_table[variant[i]].increment_n_deletion(count)
                    else:
                        self.error_table[variant[i]].increment_n_long_deletion(count)
                else: # if the long_deletion flag is on there is no updating.
                    pass
                # Set long flag deletion to True.
                long_deletion = True
            else:
                long_deletion = False

            if s == MISMATCH:
                self.error_table[variant[i]].increment_missmatch(count)

    def get_cell_text(self, library_reads):
        cell_text = []
        for letter, letter_error in self.error_table.items():
            # MissMatch, Deletion, Insertion, LongDeletion, LongInsertion
            mm, d, i, ld, li = letter_error.get_data()

            # Get letter_frequence and reads count
            l_frequency = library_reads.get_base_count_in_design(letter)
            r_count = library_reads.get_matched_reads_count()

            #
            table_output = [
                        self.get_error_rate(l_frequency, mm),
                        self.get_error_rate(l_frequency, d),
                        self.get_error_rate(l_frequency, i),
                        self.get_error_rate(r_count, ld),
                        self.get_error_rate(r_count, li)
                        ]

            # Convert values to integers
            table_output = list(map(int, table_output))

            cell_text.append(table_output)

        a_mm, a_d, a_i, a_ld, a_li = (0, 0, 0, 0, 0)
        # Sum up all errors.
        for letter, letter_error in self.error_table.items():
            # MissMatch, Deletion, Insertion, LongDeletion, LongInsertion
            mm, d, i, ld, li = letter_error.get_data()
            a_mm += mm
            a_d += d
            a_i += i
            a_ld += ld
            a_li += li

        # Get base frequency and number of reads
        b_frequency = library_reads.get_total_matched_base_count()
        r_count = library_reads.get_matched_reads_count()

        table_output = [
            self.get_error_rate(b_frequency, a_mm),
            self.get_error_rate(b_frequency, a_d),
            self.get_error_rate(b_frequency, a_i),
            self.get_error_rate(r_count, a_ld),
            self.get_error_rate(r_count, a_li)
        ]

        table_output = list(map(int, table_output))
        cell_text.append(table_output)
        return cell_text


    @staticmethod
    def get_error_rate(total, errors):
        if not errors:
            return 0
        return total / errors

    def generate_content(self, library_reads):
        col_names = ["MM rate\n[1/bases]",
                     "D rate(1nt)\n[1/bases]",
                     "I rate(1nt)\n[1/bases]",
                     "D of\nlength > 1\n[1/reads]",
                     "I of\nlength > 1\n[1/reads]"]

        row_names = ['A', 'C', 'G', 'T', 'Average']

        # Remove axis from the table.
        plt.axis('tight')
        plt.axis('off')

        # Plot table.
        cell_text = self.get_cell_text(library_reads)
        table = plt.table(cellText=cell_text,
                          rowLabels=row_names,
                          colLabels=col_names,
                          colWidths=[0.2 for _ in col_names],
                          loc='center',
                          cellLoc='center',
                          rowLoc='center',
                          bbox=[0.065, 0, 1, 1])

        table.scale(1, 3)
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        file_name = "temp/error_table.png"
        plt.savefig(file_name)
        plt.clf()

        content_array = [Content(Content.Type.TEXT, "Error Table "),
                         Content(Content.Type.IMAGE, file_name)]

        return content_array






