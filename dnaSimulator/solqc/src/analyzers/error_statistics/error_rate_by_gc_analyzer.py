from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from src.analyzers.error_statistics.deletion_analyzer import DeletionAnalyzer
from src.utils.content import Content
import pandas as pd

DELETION = "1"
INSERTION = "2"
MISMATCH = "3"

CIGAR_DICT = {
    "=": "0",
    "I": INSERTION,
    "D": DELETION,
    "X": MISMATCH
}




class ErrorRateByGCAnalyzer(DeletionAnalyzer):
    # Will store the number of deletions for each base.
    base_deletion = {}

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
        reads_df = library_reads.get_df_copy()
        reads_df['GC-content'] = reads_df.apply(func=GC_Content, axis=1, args=[library_design])
        des_len=library_design.get_design_sequences_average_len()
        reads_df['Error rate'] = reads_df.apply(func=error_rate, axis=1, args=[library_design, des_len])
        df_new = pd.DataFrame([reads_df.ix[idx]
                       for idx in reads_df.index
                       for _ in range(reads_df.ix[idx]['count'])]).reset_index(drop=True)
        medians = df_new.groupby(['GC-content'])['count'].median().values
        bp=df_new.boxplot(column=['Error rate'], by='GC-content', showmeans=True, showfliers=False, fontsize=16)
        ffff=df_new.groupby(['GC-content'])['count'].agg('sum')
        occ=ffff.values.tolist()
        totalreads=sum(occ)
        print(totalreads)
        newList = []
        for x in occ:
            newList.append(x/totalreads)
        print(occ)
        bp.plot()

        #bp.set_xticklabels(occ)
        #plt.locator_params(axis='x', tight=True)
        bp.set_xlabel("GC-content", fontsize=20)
        locs, labels = plt.xticks()
        plt.xticks(locs[::3], labels[::3], fontsize=16)
        #locs, labels = plt.xticks()
        plt.setp(labels, rotation=30, horizontalalignment='right', fontsize=16)
        '''
        #ax2 = bp.twiny()
        for tick in range(len(locs)):
            plt.text(locs[tick], medians[tick] + 0.03, occ[tick],
                     horizontalalignment='center', size='x-small', color='w', weight='semibold')

        #ax2.set_xticks([loc -1 for loc in locs[::3]])
        #ax2.set_xticklabels(occ[::3] ,rotation=90)
        #ax2.xaxis.set_ticks_position('bottom')
        #ax2.xaxis.set_ticks_position('top')
        '''
        '''
        nobs = ["n: " + str(i) for i in occ]

        # Add it to the plot
        pos = range(len(nobs))
        for tick,label in zip(pos,labels):
            plt.text(pos[tick], medians[tick] + 0.03, nobs[tick],
                horizontalalignment='center', size='x-small', color='w', weight='semibold')
            
        '''
        plt.tight_layout(pad=1)
        plt.title(" ")
        plt.suptitle("")
        plt.ylabel("Error rate", fontsize=20)
        plt.savefig("temp/boxplot2.png", dpi=300, bbox_inches = "tight")

        content_array = []
        headline = Content(Content.Type.TEXT, "Error rates by GC-Content")
        content_array.append(headline)
        content_array.append(Content(Content.Type.IMAGE, "temp/boxplot2.png"))


        return content_array







    def __str__(self):
        return 'Per base deletion analyzer'


    def update_deletion(self, path, variant):
        deletion_indices = self.locate_deletion_locations(path)

        if len(deletion_indices) == 0:
            return

        for i in deletion_indices:
            self.base_deletion[variant[i]][i] += 1

def get_error_rate(seq, path):
    counter = 0
    if float('-inf') < float(path) < float('inf'):
        for index, letter in enumerate(path):
            if letter != '0':
                counter += 1
    return counter/len(seq)

def get_del_rate(seq, path):
    counter = 0
    if float('-inf') < float(path) < float('inf'):
        for index, letter in enumerate(path):
            if letter == '1':
                counter += 1
    return counter/len(seq)

def get_cont(seq, let):
    counter = 0
    for index, letter in enumerate(seq):
        if (letter == let):
            counter += 1
    return counter/(len(seq))

def GC_Content(row, library_design):
    variant_id=row['variant_id']
    variant=library_design.get_variant_by_id(variant_id)
    variant_seq=variant.get_attribute("sequence")
    count =0.0
    for i, letter in enumerate(variant_seq):
        if(letter=='G' or letter =='C'):
            count+=1;
    return round(100*(count/len(variant_seq)))

def error_rate(row, library_design, des_len):
    cigar=row['cigar_path']
    variant_id=row['variant_id']
    variant=library_design.get_variant_by_id(variant_id)
    variant_seq=variant.get_attribute("sequence")
    des_len=len(variant_seq)
    counter=0.0
    if str(cigar) == '':
        print(cigar)
        return 0
    for i, let in enumerate(str(cigar)):
        if let != '0':
            if let=='1':
                if i-1>-0 and str(cigar)[i-1]=='1':
                    continue
            counter+=1
    return counter/des_len

def del_rate(row, library_design):
    cigar=row['cigar_path']
    counter=0.0
    if cigar.isnull() == True:
        print(cigar)
        return 0
    str_cig=str(cigar)
    for i, let in enumerate(str_cig):
        if let == '1':
            while i+1<str_cig.len() and str_cig[i+1] == '1':
                i=i+1
            counter+=1
    if counter == len(str_cig):
        print(cigar)
    return counter/len(str_cig)

def isNaN(num):
    return num != num