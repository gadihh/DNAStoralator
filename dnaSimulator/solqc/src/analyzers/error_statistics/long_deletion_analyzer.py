from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from src.analyzers.error_statistics.deletion_analyzer import DeletionAnalyzer
from src.utils.content import Content

DELETION = "1"
INSERTION = "2"
MISMATCH = "3"

CIGAR_DICT = {
    "=": "0",
    "I": INSERTION,
    "D": DELETION,
    "X": MISMATCH
}




class LongDeletionAnalyzer(DeletionAnalyzer):
    # Will store the number of deletions for each base.
    long_deletion = {
        1 : 0.0
    }
    def __init__(self):
        super(DeletionAnalyzer, self).__init__()
        self.error_rate_by_gc_cont=[]
        self.del_rate_by_gc_cont=[]

        self.longest_read_length = 0

        for key in self.base_deletion:
            self.base_deletion[key] = np.zeros(self.longest_read_length)

    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance:
            print("You can prefrom deletion analyzing on reads that did not go through edit distance")
            return None
        self.longest_read_length = library_design.get_design_longest_sequence_len()
        #reads_df: object = library_reads.get_df_copy()
        matched_reads = library_reads.get_matched_reads()
        counter_deletions=0
        read_counter =0
        total_bases=0
        for read in matched_reads:
            variant_id=read.get_variant_id()
            variant=library_design.get_variant_by_id(variant_id)
            variant_seq=variant.get_attribute("sequence")
            des_len=len(variant_seq)
            read_count = read.get_row_count()
            total_bases+=read_count*des_len
            read_counter+=read_count
            query_target_path = read.get_cigar_path()
            str_cigar = str(query_target_path)
            length = 0
            for idx, let in enumerate(str_cigar):
                if let == '1':
                    counter_deletions+=1
                    length+=1
                else:
                    if length > 0 :
                        if length in self.long_deletion:
                            self.long_deletion[length] += read_count
                        else:
                            self.long_deletion[length] = read_count
                    length=0

        for key,value in self.long_deletion.items():
            self.long_deletion[key] = value/(total_bases)

        lists = sorted(self.long_deletion.items())  # sorted by key, return a list of tuples

        x, y = zip(*lists)  # unpack a list of pairs into two tuples
        plt.clf()
        plt.plot(x, y, marker='o')
        plt.yscale("log")
        plt.grid()
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        plt.ylabel("Deletion rate", fontsize=18)
        plt.xlabel("Deletion length", fontsize=18)
        plt.savefig("temp/longDel.png", dpi=300, bbox_inches = "tight")
        #plt.show()

        content_array = []
        headline = Content(Content.Type.TEXT, "Deletion length distribution")
        content_array.append(headline)
        content_array.append(Content(Content.Type.IMAGE, "temp/longDel.png"))


        return content_array







    def __str__(self):
        return 'Per base deletion analyzer'


    def update_deletion(self, path, variant):
        deletion_indices = self.locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            self.base_deletion[variant[i]][i] += 1

