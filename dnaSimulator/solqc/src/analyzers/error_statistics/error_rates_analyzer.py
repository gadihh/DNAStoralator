from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from src.utils.content import Content
from src.analyzers.analyzer import Analyzer


DELETION = "1"
INSERTION = "2"
SUBSTITUTION = "3"





class ErrorRatesAnalyzer(Analyzer):
    # Will store the number of deletions for each base.
    base_deletion = 0.0
    long_deletion = 0.0
    substitution = 0.0
    insertion = 0.0
    total_bases = 0.0
    total_deletion = 0.0
    def __init__(self):
        pass


    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance:
            print("You can prefrom deletion analyzing on reads that did not go through edit distance")
            return None
        des_len=library_design.get_design_sequences_average_len()
        # Get all the reads we managed to match.
        matched_reads = library_reads.get_matched_reads()
        print("Number of matched reads = {}".format(len(matched_reads)))
        bar = self.get_progress_bar(len(matched_reads))
        max_error_counter = 0
        for read in matched_reads:
            variant_id=read.get_variant_id()
            variant=library_design.get_variant_by_id(variant_id)
            variant_seq=variant.get_attribute("sequence")
            des_len=len(variant_seq)
            read_count = read.get_row_count()
            query_target_path = read.get_cigar_path()
            str_cigar = str(query_target_path)
            for idx, let in enumerate(str_cigar):
                if let == '1':
                    if idx +1 < len(str_cigar) and str_cigar[idx+1]!='1':
                        if not(idx-1>=0 and str_cigar[idx-1]=='1'):
                            self.base_deletion+=read_count
                            self.total_deletion+=read_count
                    elif idx +1 < len(str_cigar) and str_cigar[idx+1]=='1':
                        if idx-1>-1 and str_cigar[idx-1] == "1":
                            self.total_deletion += read_count
                        elif idx==0 or (idx-1>-1 and str_cigar[idx-1]!="1"):
                            self.long_deletion+=read_count
                            self.total_deletion += read_count
                    else:
                        if idx-1>-1 and str_cigar[idx-1]!="1":
                            self.base_deletion += read_count
                            self.total_deletion += read_count
                    if  idx-1>-1 and str_cigar[idx-1]=="1":
                            self.total_deletion += read_count
                elif let == '2':
                    self.insertion += read_count
                elif let == '3':
                    self.substitution += read_count
            self.total_bases+=(read_count*des_len)
        print(self.base_deletion)
        self.substitution /= self.total_bases
        self.insertion /= self.total_bases
        self.base_deletion /= self.total_bases
        self.long_deletion /= self.total_bases
        self.total_deletion /= self.total_bases


        return self.generate_content()

    def __str__(self):
        return 'Per base deletion analyzer'

    def generate_content(self):
        content_array = []
        image_one_name = self.plot_error_probabilities(True)
        # Create and add content
        headline = Content(Content.Type.TEXT, "Error rates analyzer")
        content_array.append(headline)
        image_one = Content(Content.Type.IMAGE, image_one_name)
        content_array.append(image_one)

        return content_array

    def plot_error_probabilities(self, plot=True):
        labels = ['In', 'Sub', 'SingleDel', 'LongDel', 'TotalDel']
        probs = np.zeros(5).astype(float)
        probs[1] = max(self.insertion,0.0000001)
        probs[0] = max(self.substitution,0.0000001)
        probs[3] = max(self.long_deletion,0.0000001)
        probs[2] = max(self.base_deletion,0.0000001)
        probs[4] = max(self.total_deletion,0.000001)

        bar_width = 0.40
        fig = plt.figure(figsize=(8, 4))
        index = np.arange(len(labels))
        plt.yscale('log')

        color = ['#0C71B2', '#059F73', '#f08080', '#8b0000', '#FF0000']

        plt.bar(np.arange(5), probs, color=color, alpha=1.0)
        for a, b in zip(np.arange(5), probs):
            plt.text(a, b, str('%.2E' % b),  ha='center', fontsize=16)
        plt.xticks(np.arange(5), ('Sub.', 'Ins.', '1-Base Del.', 'Long Del.',  'Del.'), fontsize=16)
        plt.xlabel("Error type", fontsize=22)
        plt.ylabel("Error rate", fontsize=22)
        #plt.ylim(pow(10, -4), 1)
        plt.yticks(fontsize=16)

        image_one_name='temp/errorRates.png'
        plt.savefig(image_one_name, dpi=300, bbox_inches = "tight")


        return image_one_name


