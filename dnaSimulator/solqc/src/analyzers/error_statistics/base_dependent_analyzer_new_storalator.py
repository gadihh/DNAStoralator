from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from src.utils.content import Content
from src.analyzers.analyzer import Analyzer
import math
import json

DELETION = "1"
INSERTION = "2"
SUBSTITUTION = "3"

A=0
C=1
G=2
T=3


def BaseToInt(base):
    if base == 'A':
        return 0
    if base == 'C':
        return 1
    if base == 'G':
        return 2
    if base == 'T':
        return 3
    else:
        return 3
        #print("ERROR UNDETECTED Base")
    # exit(1)

class NewBaseDependentAnalyzer_Storalator(Analyzer):
    # Will store the number of deletions for each base.
    Del = [0.0, 0.0, 0.0, 0.0]
    Ins = [0.0, 0.0, 0.0, 0.0]
    Sub = [0.0, 0.0, 0.0, 0.0]
    LongDel = [0.0, 0.0, 0.0, 0.0]
    Ins_Bef = [0.0, 0.0, 0.0, 0.0]
    total_bases = 0.0
    total_base = [0.0, 0.0, 0.0, 0.0]
    total_long = 0
    total_sub = 0
    total_del = 0
    total_ins = 0
    total_A = 0.0
    total_C = 0.0
    total_G = 0.0
    total_T = 0.0
    def __init__(self):
        pass


    def analyze(self, library_reads, library_design):
        if not library_reads.did_edit_distance:
            print("You can prefrom deletion analyzing on reads that did not go through edit distance")
            return None

        # Get all the reads we managed to match.
        matched_reads = library_reads.get_matched_reads()
        print("Number of matched reads = {}".format(len(matched_reads)))
        bar = self.get_progress_bar(len(matched_reads))
        max_error_counter = 0
        for read in matched_reads:
            read_count = read.get_attribute("count")
            read_count = int(read_count)
            var_id = read.get_attribute("variant_id")
            read_seq = library_design.get_variant_sequence(var_id)
            for idx, let in enumerate(read_seq):
                if let=='A':
                    self.total_base[0]+=read_count
                if let=='C':
                    self.total_base[1]+=read_count
                if let=='G':
                    self.total_base[2]+=read_count
                if let=='T':
                    self.total_base[3]+=read_count
            query_target_path = read.get_cigar_path()
            variant_it = read.get_variant_id()
            variant_seq = library_design.get_variant_sequence(variant_it)
            #print(str(variant_seq))
            str_variant = str(variant_seq)
            str_cigar = str(query_target_path)
            idx_var = 0
            for idx, let in enumerate(str_cigar):
                if let == DELETION:
                    if idx +1 < len(str_cigar) and str_cigar[idx+1]!=DELETION:
                        #case of 1 base deletion
                        self.total_del+=read_count
                        if(idx_var<len(str_variant)):
                            base = BaseToInt(str_variant[idx_var])
                            self.Del[base]+=read_count
                            if (idx_var < len(read_seq)):

                                base_read = BaseToInt(str(read_seq)[idx_var])
                                #self.total_base[base_read]+=1
                                idx_var+=1
                    else:
                        #case of long deletion
                        self.total_long+=1
                        #if(idx_var<len(str_variant)):
                            #base = BaseToInt(str_variant[idx_var])
                            #self.LongDel[base]+=1
                        #idx_var+=1
                        while idx + 1 < len(str_cigar) and str_cigar[idx + 1] == '1':
                            if str_cigar[idx + 1] == '1':
                                idx += 1
                                if idx_var < len(str(read_seq)):
                                    base_read = BaseToInt(str(read_seq)[idx_var])
                                    #self.total_base[base_read] += 1
                                    idx_var+=1
                                #else:
                                #    print("here2")
                                #    base = BaseToInt(str_variant[idx_var])
                                #    self.LongDel[base]+=1
                            else:
                                base = BaseToInt(str_variant[idx_var])
                                self.LongDel[base]+=read_count
                                break
                        if idx_var < len(str(read_seq)):
                            base = BaseToInt(str_variant[idx_var])
                            self.LongDel[base]+=read_count
                elif let == INSERTION: #the ins and ins before are flipped since we have files from 5 to 3
                    self.total_ins+=read_count
                    if idx_var < len(read_seq):
                        base_read = BaseToInt(str(read_seq)[idx_var])
                        #self.total_base[base_read] += 1
                        base = BaseToInt(str(read_seq)[idx_var])
                        self.Ins_Bef[base] += read_count
                        if idx_var-1 < len(str(read_seq)) :
                            base_next = BaseToInt(str(read_seq)[idx_var-1])
                            self.Ins[base_next]+=read_count
                elif let == SUBSTITUTION:
                    self.total_sub+=read_count
                    if idx_var < len(str_variant) and idx_var < len(read_seq) :
                        base_read = BaseToInt(str(read_seq)[idx_var])
                        #self.total_base[base_read] += 1
                        base = BaseToInt(str_variant[idx_var])
                        self.Sub[base] += read_count
                    idx_var+=1

            self.total_bases+=len(str(read))
        #print(100*self.total_base[0]/self.total_bases)
        #print(100*self.total_base[1]/self.total_bases)
        #print(100*self.total_base[2]/self.total_bases)
        #print(100*self.total_base[3]/self.total_bases)


        for i in range (0,4):
            self.Del[i]/=self.total_base[i]
            #self.Del[i]*=100

            self.Sub[i]/=self.total_base[i]
            #self.Sub[i]*=100

            self.LongDel[i]/=self.total_base[i]
            #self.LongDel[i]*=100

            self.Ins[i]/=self.total_base[i]
            #self.Ins[i]*=100

            self.Ins_Bef[i] /= self.total_base[i]
            #self.Ins_Bef[i]*=100
        dic_del = {"a": self.Del[0], "c": self.Del[1], "g": self.Del[2], "t": self.Del[3]}
        dic_sub = {"a": self.Sub[0], "c": self.Sub[1], "g": self.Sub[2], "t": self.Sub[3]}
        dic_long = {"a": self.LongDel[0], "c": self.LongDel[1], "g": self.LongDel[2], "t": self.LongDel[3]}
        dic_ins = {"a": self.Ins[0], "c": self.Ins[1], "g": self.Ins[2], "t": self.Ins[3]}
        dic_pre_ins = {"a": self.Ins_Bef[0], "c": self.Ins_Bef[1], "g": self.Ins_Bef[2], "t": self.Ins_Bef[3]}
        dic_output = {
            "one_base_delete" : dic_del,
            "substitution" : dic_sub,
            "long_delete" : dic_long,
            "insert" : dic_ins,
            "pre_insert": dic_pre_ins
        }
        with open("temp/storalator_config.json") as f:
            data = json.load(f)
        data.update(dic_output)
        with open("temp/storalator_config.json", 'w') as f:
            json.dump(data, f, indent=4, separators=(',', ': '))
        return self.generate_content()

    def __str__(self):
        return 'Per base deletion analyzer'

    def generate_content(self):
        content_array = []

        image_one_name = self.heat_map_plot()
        # Create and add content
        image_one = Content(Content.Type.IMAGE, image_one_name)
        content_array.append(image_one)

        return content_array


    def heat_map_plot(self):
        bases = ["A", "C", "G", "T"]
        errors = ["Substitution", "Inserted Symbol", "Symbol Pre-Insertion", "1-Base Deletion", "Long Deletion"]
        rates = np.array([[self.Sub[0], self.Sub[1], self.Sub[2], self.Sub[3]],
                          [self.Ins[0], self.Ins[1], self.Ins[2], self.Ins[3]],
                          [self.Ins_Bef[0], self.Ins_Bef[1], self.Ins_Bef[2], self.Ins_Bef[3]],
                          [self.Del[0], self.Del[1], self.Del[2], self.Del[3]],
                          [self.LongDel[0], self.LongDel[1], self.LongDel[2], self.LongDel[3]]])

        fig, ax = plt.subplots()
        im = ax.imshow(rates)
        # We want to show all ticks...
        ax.set_xticks(np.arange(len(bases)))
        ax.set_yticks(np.arange(len(errors)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(bases, fontsize=28)
        ax.set_yticklabels(errors, fontsize=28)
        colors=["k", "w", "w", "k", "w"]
        for i in range(len(bases)):
            for j in range(len(errors)):
                ext = ax.text(i, j, round(rates[j, i], 4),
                              ha="center", va="center", color=colors[j], fontsize=20)
        cbar = ax.figure.colorbar(im, ax=ax, **{})
        cbar.ax.set_ylabel("Error rates presented in percents", rotation=-90, va="bottom", fontsize=28)
        #plt.title("Symbol Dependent Errors")
        image_one_name='temp/SymbolDependent.png'
        plt.savefig(image_one_name)
        plt.tight_layout(pad=0.8)
        #plt.show()

        return image_one_name



