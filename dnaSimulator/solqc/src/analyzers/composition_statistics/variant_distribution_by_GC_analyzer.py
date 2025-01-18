import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import seaborn as sns
import pandas as pd
import math as math

import pickle # TODO : remove after tali haran data as been processed.

from src.analyzers.analyzer import Analyzer
from src.utils.content import Content


class VariantDistributionByGCAnalyzer(Analyzer):
    name = "variant distribution analyzer"
    content_array = []
    stat = []
    title = []
    stat_sorted = []
    title_sorted = []
    def __init__(self):
        pass

    def analyze(self, library_reads, library_design):
        if not library_reads.did_matching():
            print("You can not preform variant distribution analysis on a library that did not do matching.")
            return
        # Get library reads dataframe.
        plt.clf()

        reads_df = library_reads.get_df_copy()
        reads_df['GC_cont'] = reads_df.apply(func=GC_Content, axis=1, args=[library_design])
        print(reads_df['GC_cont'].max())
        print(reads_df['GC_cont'].min())

        max_GC=math.ceil(reads_df['GC_cont'].max())
        min_GC=int(reads_df['GC_cont'].min())
        jump=((max_GC-min_GC)/5)
        # Get variant id series and remove the rows where variant id is unknown.
        group_by_variant = reads_df.groupby(pd.cut(reads_df['GC_cont'], np.arange(min_GC, max_GC+jump, jump)))
        last_GC_cont = []

        for GC_con, group in group_by_variant:
            if(len(group)>0):

                # Count the number of reads each variant has.
                variants_count=group.groupby(['variant_id'])['count'].agg('sum')
                self.generate_content({'variants_count': variants_count}, GC_con, self.content_array, self.stat, self.title)
                last_GC_cont=GC_con
        self.content_array.append(Content(Content.Type.TEXT, "Sorted bar plot of the number of reads per variant, stratified by GC content"))
        self.content_array.append(Content(Content.Type.IMAGE, "temp/variant_distribution{}.png".format(last_GC_cont)))
        plt.clf()
        
        i, j, d = (0,0,0)
        x = math.ceil(math.sqrt(len(self.stat)))
        fig, ax = plt.subplots(x, x, sharey=True,sharex=True)
        for idx, sorted_data in enumerate(self.stat):
            ax[i][j].bar(np.arange(len(sorted_data)), sorted_data, alpha=1) # Edge color solution is a bit hacky but I don't have a better idea
            ax[i][j].set_title("{}".format(self.title[idx]), fontsize=11)
            if d % 2 == 0:
                i=(i+1)%x
                if i == 0:
                    d+=1
                    j = (j + 1) % x
                    d += 1
            else:
                j=(j+1)%x
                i=0
                d+=1
        plt.tight_layout()
        plt.savefig("variant_distributionTotal.png")
        self.content_array.append(Content(Content.Type.IMAGE, "variant_distributionTotal.png"))

        plt.clf()

        for GC_con, group in group_by_variant:
            if (len(group) > 0):
                # Count the number of reads each variant has.
                variants_count=group.groupby(['variant_id'])['count'].agg('sum')

                self.generate_content_sorted({'variants_count': variants_count}, GC_con, self.content_array, self.stat_sorted,
                                             self.title_sorted)
                last_GC_cont = GC_con
        self.content_array.append(Content(Content.Type.NEW_PAGE, ""))
        self.content_array.append(Content(Content.Type.TEXT, "Histogram of the number of reads, stratified by GC Content"))
        self.content_array.append(Content(Content.Type.IMAGE, "temp/variant_distribution_sorted{}.png".format(last_GC_cont)))
        plt.clf()
        i, j, d = (0, 0, 0)
        x = math.ceil(math.sqrt(len(self.stat_sorted)))
        fig, ax = plt.subplots(x, x, sharey=True, sharex=True)
        for idx, sorted_data in enumerate(self.stat_sorted):
            ax[i][j].bar(np.arange(len(sorted_data)), sorted_data,
                         alpha=1)  # Edge color solution is a bit hacky but I don't have a better idea
            ax[i][j].set_title("{}".format(self.title_sorted[idx]), fontsize=11)
            if d % 2 == 0:
                i = (i + 1) % x
                if i == 0:
                    d += 1
                    j = (j + 1) % x
                    d += 1
            else:
                j = (j + 1) % x
                i = 0
                d += 1
        plt.savefig("temp/variant_distribution_sorted_total.png")
        self.content_array.append(Content(Content.Type.IMAGE, "temp/variant_distribution_sorted_total.png"))
        plt.clf()
        # save data.
        return self.content_array

    @staticmethod
    def generate_content(data, GC_Con, content_array, stat, title):
        variants_count = data['variants_count']

        var_distribution_data = variants_count

        # Generating the number of reads per variant graph.
        plt.xlabel("Variant rank", fontsize=20)
        plt.ylabel("Number of reads", fontsize=20)
        sorted_data = np.flipud(np.sort(var_distribution_data.values))
        label = "{}".format(GC_Con)
        plt.bar(np.arange(len(sorted_data)), sorted_data, alpha=0.4, label=label, align='center')
        plt.plot(np.arange(len(sorted_data)), sorted_data, linewidth=0.5)

        plt.grid()
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16, rotation=30)
        plt.legend(fontsize=16)

        plt.savefig("temp/variant_distribution{}.png".format(GC_Con), dpi=300, bbox_inches = "tight")
        image_one = Content(Content.Type.TEXT, "Variants Distribution for variants with GC Content {}".format(GC_Con))
        stat.append(sorted_data)
        title.append(label)
        return content_array

    @staticmethod
    def generate_content_sorted(data, GC_Con, content_array, stats_sorted, title_sorted):
        variants_count = data['variants_count']

        var_distribution_data = variants_count


        # Generating the number of reads distribution.
        plt.xlabel("Number of reads", fontsize=22)
        plt.ylabel("Number of variants", fontsize=22)

        # Set the number of bins to the number of distinct values.
        bins = len(set(variants_count.values))

        # Set the y-axis range. (otherwise ymax is the len of the data and the values ares almost invisible).
        most_common, num_most_common = Counter(variants_count.values).most_common(1)[0]
        print("Most common value is {} and it appeared {} times".format(most_common, num_most_common))
        label = "{}".format(GC_Con)
        labels, counts = np.unique(variants_count.values, return_counts=True)
        plt.bar(labels, counts, align='center', alpha=0.4, label=label)
        plt.plot(labels, counts, linewidth=0.5)

        plt.legend(fontsize=18)
        plt.grid()
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)

        plt.savefig("temp/variant_distribution_sorted{}.png".format(GC_Con), dpi=300, bbox_inches = "tight")
        stats_sorted.append(counts)
        title_sorted.append(label)

        return content_array

def GC_Content(row, library_design):
    variant_id=row['variant_id']
    variant=library_design.get_variant_by_id(variant_id)
    variant_seq=variant.get_attribute("sequence")
    count =0.0
    for i, letter in enumerate(variant_seq):
        if(letter=='G' or letter =='C'):
            count+=1;
    return 100*(count/len(variant_seq))