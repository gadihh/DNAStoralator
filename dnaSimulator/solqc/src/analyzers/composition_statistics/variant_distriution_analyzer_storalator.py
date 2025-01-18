'''
Variant Distribution Analyzer
The additional section will include a grphical representation of the distribution
of the number of reads per variant, as well as the number of reads each variant has.
'''

import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import seaborn as sns

from src.analyzers.analyzer import Analyzer
from src.utils.content import Content
import json


class VariantDistributionAnalyzer_Storalator(Analyzer):
    name = "variant distribution analyzer storalator"
    display_rank = 100

    def __init__(self):
        pass

    def analyze(self, library_reads, library_design):
        if not library_reads.did_matching():
            print("You can not preform variant distribution analysis on a library that did not do matching.")
            return

        # Get library reads dataframe.
        reads_df = library_reads.get_matched_dataframe()
        # Get variant id series and remove the rows where variant id is unknown.
        variants_count = reads_df.groupby(['variant_id'])['count'].agg('sum').values
        # # Count the number of reads each variant has.

        # save data.
        return self.generate_content({'variants_count': variants_count})

    @staticmethod
    def generate_content(data):
        variants_count = data['variants_count']
        var_distribution_data = variants_count
        plt.clf()
        # Generating the number of reads per variant graph.
        plt.xlabel("Variant rank", fontsize=20)
        plt.ylabel("Number of reads", fontsize=20)
        plt.xticks(fontsize=16, rotation=30)
        plt.yticks(fontsize=16)
        plt.ylim(0, np.max(var_distribution_data) + 2)
        sorted_data = np.flipud(np.sort(var_distribution_data))
        plt.bar(np.arange(len(sorted_data)), sorted_data, edgecolor=sns.color_palette()[0]) # Edge color solution is a bit hacky but I don't have a better idea.
        plt.grid()
        plt.savefig("temp/variant_distribution.png", dpi=300, bbox_inches = "tight")

        #print(list(x_axis_val))
        plt.clf()

        # Adding the unmatched reads

        # Generating the number of reads distribution.
        plt.xlabel("Number of reads", fontsize=20)
        plt.ylabel("Number of variants", fontsize=20)

        # Set the number of bins to the number of distinct values.
        bins = len(set(variants_count))

        # Set the y-axis range. (otherwise ymax is the len of the data and the values ares almost invisible).
        most_common, num_most_common = Counter(variants_count).most_common(1)[0]
        print("Most common value is {} and it appeared {} times".format(most_common, num_most_common))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        labels, counts = np.unique(variants_count, return_counts=True)
        plt.bar(labels, counts, align='center')
        plt.grid()
        plt.savefig("temp/variant_read_distribution.png", dpi=300, bbox_inches = "tight")
        x_val_explicit_vector = list((np.arange(len(sorted_data))))

        for group in [x_val_explicit_vector, sorted_data, labels, counts]:    
            for index, item in enumerate(group):
                group[index] = int(item)
        
        dic_output = {
            "x_val_explicit_vector" : x_val_explicit_vector,
            "y_val_explicit_vector" : sorted_data.tolist(),
            "x_val_histo" : labels.tolist(),
            "y_val_histo" : counts.tolist()
        }

        f = open("temp/storalator_config_cluster_size.json", "w")
        json.dump(dic_output, f, indent = 6)
        f.close()
        plt.clf()

        content = [Content(Content.Type.TEXT, "Sorted bar plot of the number of filtered reads per variant"),
                   Content(Content.Type.IMAGE, "temp/variant_distribution.png"),
                   Content(Content.Type.NEW_PAGE, ""),
                   Content(Content.Type.TEXT, "Histogram of the number of filtered reads per variant"),
                   Content(Content.Type.IMAGE, "temp/variant_read_distribution.png")]

        return content

