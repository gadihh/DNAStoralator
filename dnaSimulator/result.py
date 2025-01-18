import os
import re
from PyQt5.QtWidgets import *
import reconResult_ui
from PyQt5.QtGui import QPixmap
import json


class ReconResult(QMainWindow, reconResult_ui.Ui_recon_result_window):
    def __init__(self, json_results):
        QMainWindow.__init__(self)
        self.setupUi(self)
        input_file = json_results.__dict__['input_file']
        self.setWindowTitle(f'Reconstruction Result for input: {input_file}')
        self.json_results = json_results
        self.result_img_path = json_results.__dict__['hist_file']
        for json_field, recon_field in zip(json_result_to_recon_result.keys(), json_result_to_recon_result.values()):
            self.__dict__[recon_field].setText(json_results.__dict__[json_field])
        self.insert_img(self.recon_histogram_img, self.result_img_path)

    def open_window(self):
        self.show()

    def insert_img(self, img_label, img_path):
        img_label.setPixmap(QPixmap(img_path))
        img_label.show()


json_result_to_recon_result = {'error_rate': 'recon_error_rate_val',
                               'substitution_rate': 'recon_substitution_rate_val',
                               'deletion_rate': 'recon_deletion_rate_val',
                               'insertion_rate': 'recon_insertion_rate_val',
                               'success_rate': 'recon_success_rate_val',
                               'num_of_clusters': 'recon_clusters_val',
                               'maj_test': 'recon_maj_test_val'}


class JsonReconResult:
    def __init__(self):
        self.input_file = ''
        self.output_file = ''
        self.algo = ''
        self.error_rate = ''
        self.deletion_rate = ''
        self.insertion_rate = ''
        self.substitution_rate = ''
        self.success_rate = ''
        self.num_of_clusters = ''
        self.maj_test = ''
        self.hist_file = ''


class ClusteringResult:
    def __init__(self, name, index, tech, seconds, false_positives, false_negatives, thrown_strands):
        self.name = name
        self.index = index
        self.tech = tech
        self.seconds = seconds
        self.false_positives = false_positives
        self.false_negatives = false_negatives
        self.thrown_strands = thrown_strands


def parse_recon_output(hist_folder, input_file_name, output_file_name, algo, raw_output):
    parse_dict = {
        'Substitution': 'substitution_rate',
        'Deletion': 'deletion_rate',
        'Insertion': 'insertion_rate',
        'Error': 'error_rate',
        'Success': 'success_rate',
        'clusters': 'num_of_clusters',
        'test': 'maj_test'
    }

    json_recon_result = JsonReconResult()

    json_recon_result.__dict__['input_file'] = input_file_name
    json_recon_result.__dict__['output_file'] = output_file_name
    json_recon_result.__dict__['algo'] = algo

    keywords_in_result = [key for key in parse_dict.keys() if key in raw_output]
    for key_word in keywords_in_result:
        key_split = re.split(key_word, raw_output)[1]
        found = re.search(r"[-+]?(?:\d*\.\d+|\d+)", key_split)
        field_value = found.group(0) if found is not None else "nan"
        json_recon_result.__dict__[parse_dict[key_word]] = field_value
    for key_word in filter(lambda key: key not in keywords_in_result, parse_dict.keys()):
        json_recon_result.__dict__[parse_dict[key_word]] = "nan"
    json_recon_result.hist_file = create_recon_result_histogram(hist_folder, input_file_name, output_file_name)
    assert (json_recon_result.hist_file != "")
    return json_recon_result


def parse_hist_results(hist_folder, recon_result_path, recon_result_file_name):
    print("output:" + recon_result_path)
    num_clusters = 0
    start_copying = 0
    x = []
    y = []
    hist_aux_path = hist_folder + "/" + re.split(".txt", recon_result_file_name)[0] + "_hist.txt"
    source = open(recon_result_path, 'r')
    f = open(hist_aux_path, 'w', newline='\n')

    for line in source:
        if line.find('rate') > 0:
            start_copying = 0
        if line.find('hist') > 0:
            start_copying = 1
        if line.find('clusters') > 0:
            numbers = [int(s) for s in re.findall(r'\b\d+\b', line)]
            num_clusters = numbers[0]
        if start_copying == 1 and not line.find('hist') > 0:
            f.write(line)

    source.close()
    f.close()

    source = open(hist_aux_path, 'r')
    for line in source:
        if (line == "\n"):
            continue
        line_list = re.split(r'\t+|\s+', line)  # seperates a string to a list by a delimeter of spaces or tabs
        index = int(line_list[0].strip())
        value = int(line_list[1].strip())
        x.insert(index, index)
        y.insert(index, (value / num_clusters) * 100)
        # y[index] = (value / num_clusters) * 100
    source.close()
    os.remove(hist_aux_path)
    return x, y, num_clusters


def create_recon_result_histogram(hist_folder, recon_input_path, recon_result_path):
    import matplotlib.pyplot as plt
    recon_result_file_name = recon_result_path.split("/")[-1]
    hist_path = hist_folder + "/" + re.split(".txt", recon_result_file_name)[0] + "_hist.png"
    x, y, num_clusters = parse_hist_results(hist_folder, recon_result_path, recon_result_file_name)

    plt.xticks(x)

    plt.figure()  # creates a new plot, so it doesn't plot a couple in one figure
    plt.scatter(x, y, color='r', zorder=2)
    plt.plot(x, y, color='b', zorder=1)

    plt.title("Reconstruction Errors Histogram")
    plt.xlabel("Number of edit errors")
    plt.ylabel("Fraction of reads")

    plt.savefig(hist_path)
    return hist_path


def update_json_file(json_path, res_dict, search_keys):
    """
    Update json entry specified in keys dict with a res_dict values. If such entry does not exist,
    a new one will be appended.
    Passing res_dict as None equals to deleting all entries that fit the search keys.
    """
    with open(json_path) as rf:
        new_entries = []
        old_entries = json.load(rf)
        # filter out irrelevant entries according to the search_keys:
        for entry in old_entries:
            to_add = False
            for key in search_keys:
                to_add = to_add + bool(entry[key] != search_keys[key])
            if to_add:
                new_entries.append(entry)
        # add the relevant entry:
        if res_dict is not None:
            new_entries.append(res_dict)
    json_dump(json_path, new_entries)


def rename_json_entry(json_path, old_name, new_name, search_key):
    with open(json_path) as rf:
        entries = json.load(rf)
        for entry in entries:
            if entry[search_key] == old_name:
                entry[search_key] = new_name
    json_dump(json_path, entries)


def json_dump(path, entries):
    with open(path, 'w') as wf:
        json.dump(entries, wf, indent=4, separators=(',', ': '))
