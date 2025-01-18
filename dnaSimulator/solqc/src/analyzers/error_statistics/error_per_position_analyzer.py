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

class ErrorPerPositionAnalyzer(DeletionAnalyzer):
    # Will store the number of deletions for each base.
    single_deletion = {
        1 : 0.0
    }
    insertion = {
        1 : 0.0
    }
    substitution = {
        1 : 0.0
    }
    long_deletion = {
        1 : 0.0
    }
    total_deletion = {
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
        matched_reads = library_reads.get_matched_reads()
        counter_reads=0
        des_len=library_design.get_design_longest_sequence_len()
        for read in matched_reads:
            read_count = read.get_attribute("count")
            counter_reads+=read_count
            query_target_path = read.get_cigar_path()
            str_cigar = str(query_target_path)
            for idx, let in enumerate(str_cigar):
                length = 0
                if let == '1':
                    length+=1
                    if idx + 1 < len(str_cigar) and str_cigar[idx + 1] != "1":
                        if not(idx-1>0 and str_cigar[idx-1]=="1"): #only if we are not in the end of long deletion
                            if idx in self.single_deletion:
                                self.single_deletion[idx] +=read_count
                            else:
                                self.single_deletion[idx] =read_count
                            if idx in self.total_deletion:
                                self.total_deletion[idx] +=read_count
                            else:
                                self.total_deletion[idx] =read_count
                    elif idx +1 <len(str_cigar) and str_cigar[idx + 1] == "1":
                        if idx-1>-1 and str_cigar[idx-1] == "1":
                            if idx in self.total_deletion:
                                self.total_deletion[idx] += read_count
                            else:
                                self.total_deletion[idx] = read_count
                        elif idx==0 or (idx-1>-1 and str_cigar[idx-1]!="1"):
                            if idx in self.long_deletion:
                                self.long_deletion[idx] +=read_count
                            else:
                                self.long_deletion[idx] =read_count
                            if idx in self.total_deletion:
                                self.total_deletion[idx] += read_count
                            else:
                                self.total_deletion[idx] = read_count
                if let == '2':
                    if idx in self.insertion:
                        self.insertion[idx]+=read_count
                    else:
                        self.insertion[idx]=read_count
                if let == '3':
                    if idx in self.substitution:
                        self.substitution[idx]+=read_count
                    else:
                        self.substitution[idx]=read_count

        for key, value in self.long_deletion.items():
            self.long_deletion[key] = value/counter_reads
        for key, value in self.single_deletion.items():
            self.single_deletion[key] = value/counter_reads
        for key, value in self.insertion.items():
                self.insertion[key] = value / counter_reads
        for key, value in self.substitution.items():
                self.substitution[key] = value / counter_reads
        for key, value in self.total_deletion.items():
            if value > counter_reads:
                print(value)
                print(key)
            self.total_deletion[key] = value / counter_reads
        lists = sorted(self.substitution.items())  # sorted by key, return a list of tuples
        plt.clf()
        x, y = zip(*lists)  # unpack a list of pairs into two tuples
        plt.plot(x[0:des_len], y[0:des_len], color='b', label="Sub.")

        lists = sorted(self.insertion.items())  # sorted by key, return a list of tuples

        x, y = zip(*lists)  # unpack a list of pairs into two tuples

        plt.plot(x[0:des_len], y[0:des_len], color='g', label="Ins.")

        lists = sorted(self.single_deletion.items())  # sorted by key, return a list of tuples

        x, y = zip(*lists)  # unpack a list of pairs into two tuples

        plt.plot(x[0:des_len], y[0:des_len], color='r', label="1-Base Del.")

        lists = sorted(self.long_deletion.items())  # sorted by key, return a list of tuples

        x, y = zip(*lists)  # unpack a list of pairs into two tuples

        plt.plot(x[0:des_len], y[0:des_len], color='orange', label="Long Del.")

        lists = sorted(self.total_deletion.items())  # sorted by key, return a list of tuples

        x, y = zip(*lists)  # unpack a list of pairs into two tuples

        plt.plot(x[0:des_len], y[0:des_len], color='black', label="Total Del.")
        plt.xlim(0, des_len-1)
        print(des_len)
        plt.grid()
        plt.yscale('log')
        plt.legend(fontsize=12)
        plt.xlabel("Position", fontsize=20)
        plt.ylabel("Error rate", fontsize=20)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        #plt.tight_layout()
        plt.savefig("temp/PerPosition.png", dpi=300, bbox_inches = "tight")
        #plt.show()

        content_array = []
        headline = Content(Content.Type.TEXT, "Error rates per position (5' to 3')")
        content_array.append(headline)
        content_array.append(Content(Content.Type.IMAGE, "temp/PerPosition.png"))


        return content_array



    def __str__(self):
        return 'Per base deletion analyzer'


    def update_deletion(self, path, variant):
        deletion_indices = self.locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            self.base_deletion[variant[i]][i] += 1

