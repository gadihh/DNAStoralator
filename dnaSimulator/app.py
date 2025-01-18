import os
import re
import shlex
import subprocess
import time
from functools import partial
import ast
import time
from datetime import datetime
from worker import *
from PyQt5 import QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import *
import sys
import json
import platform
import functools as ft
import dnaSimulator_ui2
from simulator import *
from index_clustering import *
from hash_based_clustering import *
from result import *
from filepath import *
import solqc_help_window
from table import *
from genericpath import isdir

class dnaSimulator(QMainWindow, dnaSimulator_ui2.Ui_dnaSimulator):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)

        # set the title
        self.setWindowTitle('DNA Storalator')
        self.inputDNAPath = ""

        # initialize general errors
        self.general_errors = {
            'd': '',
            'ld': '',
            'i': '',
            's': ''
        }

        # initialize per-base errors
        self.per_base_errors = {
            'A': {'s': '', 'i': '', 'pi': '', 'd': '', 'ld': ''},
            'C': {'s': '', 'i': '', 'pi': '', 'd': '', 'ld': ''},
            'G': {'s': '', 'i': '', 'pi': '', 'd': 'test', 'ld': ''},
            'T': {'s': '', 'i': '', 'pi': '', 'd': '', 'ld': ''}
        }

        self.dist_info = {
            'type': '',
            'value': '',
            'min': 0,
            'max': 0
        }
        self.chosen_tech_dict = {'Twist Bioscience + Ilumina miSeq': 'miseq_twist',
                                 'CustomArray + Ilumina miSeq': 'miseq_custom',
                                 'Twist Bioscience + Ilumina NextSeq': 'nextseq_twist',
                                 'Integrated DNA Technology (IDT) + MinION': 'minion_idt',
                                 'Stutter': 'stutter_def',
                                 'Other': 'other'}

        self.recon_algo_dict = {
            'Hybrid Reconstruction Algorithm': 'Hybrid',
            'Divider BMA Reconstruction Algorithm': 'DivBMA',
            'BMA Look Ahead Reconstruction Algorithm': 'BMALookahead',
            'Iterative Reconstruction Algorithm': 'Iterative',
            'Mock Reconstruction Algorithm': 'mockReconstruction'
        }
        # list of tables
        self.tables = [self.error_sim_table, self.clustering_table, self.reconstruction_table]
        self.progress_tables = [self.clustering_progress_table, self.recon_progress_table]

        self.error_sim_table_cols = {'name': 0, 'size': 1, 'time': 2, 'action': 3, 'progress': 4}
        self.cluster_table_cols = {'name': 0, 'size': 1, 'time': 2, 'result': 3, 'action': 4, '': 5}
        self.recon_table_cols = {'name': 0, 'size': 1, 'time': 2, 'result': 3, 'action': 4, '': 5}

        self.cluster_progress_cols = {'name': 0, 'tech': 1, 'index': 2, 'stop': 3, 'progress': 4}
        self.recon_progress_cols = {'name': 0, 'algo': 1, 'stop': 2, 'progress': 3}

        self.table_cols = [self.error_sim_table_cols, self.cluster_table_cols, self.recon_table_cols]
        self.progress_table_cols = [self.cluster_progress_cols, self.recon_progress_cols]

        self.set_working_dirs()
        self.set_json()

        self.barcode_start = 0
        self.barcode_end = 0
        self.max_clustering_edit_dist = 0

        self.minHash_local_steps = 220
        self.minHash_min_local_steps = 20
        self.minHash_max_local_steps = 1000

        self.minHash_similarity_threshold = 15
        self.minHash_min_similarity_threshold = 1
        self.minHash_max_similarity_threshold = 1000

        self.minHash_window_size = 4
        self.minHash_min_window_size = 1
        self.minHash_max_window_size = 1000

        self.reconstruction_algo = ''
        self.clustering_algo = ''
        self.error_chosen_technology = ''
        self.cluster_chosen_technology = ''
        self.recon_chosen_technology = ''
        self.clustering_index = None

        self.clustering_orig_input = ''
        self.recon_input = ''
        self.clustering_output = ''
        self.error_output = ''
        self.shuffled_output = ''

        self.clustering_tab_name = 'Clustering'
        self.reconstruction_tab_name = 'Reconstruction'

        self.cluster_workers = []
        self.error_workers = []
        self.recon_processes = []

        self.recon_result_windows = []
        self.set_minHash_cluster_values()

        # connect push buttons to an event
        self.browse_PushButton.clicked.connect(self.openFileDialog)

        self.run_error_simulator_PushButton.clicked.connect(self.runErrorSimulator)
        self.run_clustering_pushButton.clicked.connect(self.runClustering)
        self.reconstruction_run_pushButton.clicked.connect(self.run_reconstruction_algo)

        # connect radio buttons to an event
        self.Ilumina_miSeq_radioButton.toggled.connect(self.ilumina_miSeq_chosen)
        self.Ilumina_NextSeq_radioButton.toggled.connect(self.ilumina_NextSeq_chosen)
        self.MinION_radioButton.toggled.connect(self.MinIon_chosen)
        self.twist_bioscience_radioButton.toggled.connect(self.twist_bioscience_chosen)
        self.customArray_radioButton.toggled.connect(self.customArray_chosen)
        self.IDT_radioButton.toggled.connect(self.IDT_chosen)
        self.stutter_radioButton.toggled.connect(self.Stutter_chosen)
        self.default_radioButton.toggled.connect(self.Default_chosen)
        self.cont_radioButton.toggled.connect(self.continuous_amount)
        self.vector_radioButton.toggled.connect(self.vector_amount)

        self.substitution_doubleSpinBox.textChanged.connect(self.set_substitution)
        self.insertion_doubleSpinBox.textChanged.connect(self.set_insertion)
        self.one_base_del_doubleSpinBox.textChanged.connect(self.set_one_base_del)
        self.long_del_doubleSpinBox.textChanged.connect(self.set_long_del)

        # connect the spinboxes of per base errors
        self.A_substitution_doubleSpinBox.textChanged.connect(self.set_A_substitution)
        self.C_substitution_doubleSpinBox.textChanged.connect(self.set_C_substitution)
        self.G_substitution_doubleSpinBox.textChanged.connect(self.set_G_substitution)
        self.T_substitution_doubleSpinBox.textChanged.connect(self.set_T_substitution)

        self.A_insertion_doubleSpinBox.textChanged.connect(self.set_A_insertion)
        self.C_insertion_doubleSpinBox.textChanged.connect(self.set_C_insertion)
        self.G_insertion_doubleSpinBox.textChanged.connect(self.set_G_insertion)
        self.T_insertion_doubleSpinBox.textChanged.connect(self.set_T_insertion)

        self.A_pre_insertion_doubleSpinBox.textChanged.connect(self.set_A_pre_insertion)
        self.C_pre_insertion_doubleSpinBox.textChanged.connect(self.set_C_pre_insertion)
        self.G_pre_insertion_doubleSpinBox.textChanged.connect(self.set_G_pre_insertion)
        self.T_pre_insertion_doubleSpinBox.textChanged.connect(self.set_T_pre_insertion)

        self.A_one_base_del_doubleSpinBox.textChanged.connect(self.set_A_one_base_del)
        self.C_one_base_del_doubleSpinBox.textChanged.connect(self.set_C_one_base_del)
        self.G_one_base_del_doubleSpinBox.textChanged.connect(self.set_G_one_base_del)
        self.T_one_base_del_doubleSpinBox.textChanged.connect(self.set_T_one_base_del)

        self.A_long_del_doubleSpinBox.textChanged.connect(self.set_A_long_del)
        self.C_long_del_doubleSpinBox.textChanged.connect(self.set_C_long_del)
        self.G_long_del_doubleSpinBox.textChanged.connect(self.set_G_long_del)
        self.T_long_del_doubleSpinBox.textChanged.connect(self.set_T_long_del)

        self.barcode_start_spinBox.textChanged.connect(self.set_barcode_start)
        self.barcode_end_spinBox.textChanged.connect(self.set_barcode_end)
        self.max_edit_dist_spinBox.textChanged.connect(self.set_clustering_edit_dist)

        self.min_amount_spinBox.textChanged.connect(self.set_min_amount)
        self.max_amount_spinBox.textChanged.connect(self.set_max_amount)
        self.value_lineEdit.textChanged.connect(self.set_value_amount)

        self.minHash_local_steps_spinBox.textChanged.connect(self.set_minHash_local_steps)
        self.minHash_window_size_spinBox.textChanged.connect(self.set_minHash_window_size)
        self.minHash_similarity_threshold_spinBox.textChanged.connect(self.set_minHash_similarity_threshold)

        self.reconstruction_listWidget.addItem('Hybrid Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Divider BMA Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('BMA Look Ahead Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Iterative Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Mock Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Not a real algo')

        self.reconstruction_listWidget.setCurrentRow(0)
        self.reconstruction_algo = self.reconstruction_listWidget.item(0).text()
        self.reconstruction_listWidget.currentItemChanged.connect(self.set_reconstruction_algo)

      
        self.stackedWidget.hide()
        self.clustering_settings_label.setVisible(False)

        self.technology_comboBox.addItem('Twist Bioscience + Ilumina miSeq')
        self.technology_comboBox.addItem('CustomArray + Ilumina miSeq')
        self.technology_comboBox.addItem('Twist Bioscience + Ilumina NextSeq')
        self.technology_comboBox.addItem('Integrated DNA Technology (IDT) + MinION')
        self.technology_comboBox.currentTextChanged.connect(self.set_chosen_clustering_technology)

        self.minHash_technology_comboBox.addItem('Twist Bioscience + Ilumina miSeq')
        self.minHash_technology_comboBox.addItem('CustomArray + Ilumina miSeq')
        self.minHash_technology_comboBox.addItem('Twist Bioscience + Ilumina NextSeq')
        self.minHash_technology_comboBox.addItem('Integrated DNA Technology (IDT) + MinION')
        self.minHash_technology_comboBox.currentTextChanged.connect(self.set_chosen_minHash_clustering_technology)

        self.clustering_index_comboBox.addItem('4')
        self.clustering_index_comboBox.addItem('5')
        self.clustering_index_comboBox.addItem('9')
        self.clustering_index_comboBox.currentTextChanged.connect(self.set_clustering_index)

        self.clustering_listWidget.currentItemChanged.connect(self.set_clustering_algo)

        self.cont_radioButton.setVisible(False)
        self.vector_radioButton.setVisible(False)
        self.value_lineEdit.setVisible(False)
        self.min_amount_spinBox.setVisible(False)
        self.max_amount_spinBox.setVisible(False)
        self.value_label.setVisible(False)
        self.min_label.setVisible(False)
        self.max_label.setVisible(False)
        self.user_defined_copied_checkBox.stateChanged.connect(self.user_defined_amount)

        # add tool tips in UI
        self.value_lineEdit.setToolTip('Vector syntax example: [1, 10, 5, 9, 6, 200]\n\n'
                                       'Continuous function should be python syntax function. For example: x ** 2')
        self.file_path_lineEdit.setToolTip('A path to your strands file')

        # table work
        self.init_tables()
        self.fill_tables()

        # initiate clustering tab
        self.clustering_listWidget.setCurrentRow(1)
        self.set_clustering_algo()

        # *** SOLQC ***
        # solqc inputs initializations
        self.design_file_path = ""
        self.seq_file_path = ""
        self.config_file_path = ""
        self.output_folder_path = ""
        self.last_dir = "./"
        self.config_file_path_line_edit.setReadOnly(True)
        self.design_file_path_line_edit.setReadOnly(True)
        self.output_file_path_line_edit.setReadOnly(True)
        self.sequencing_file_path_line_edit.setReadOnly(True)
        self.design_file_path_button.clicked.connect(self.open_design_file_dialog)
        self.sequencing_file_path_button.clicked.connect(self.open_sequencing_file_dialog)
        self.config_file_path_button.clicked.connect(self.open_config_file_dialog)
        self.output_file_path_button.clicked.connect(self.open_output_folder_dialog)
        self.start_solqc_button.clicked.connect(self.start_solqc)
        self.json_pushBrowse.clicked.connect(self.search_for_json_file)
        self.json_pushLoad.clicked.connect(self.open_and_load_json_file)

        # solqc help window
        self.solqc_help_window = solqcHelpWindow()
        # self.solqc_help_window_visible = False
        # self.solqc_help_window.hide()
        self.solqc_help_button.clicked.connect(self.show_solqc_help)

        # solqc optional synth & seq techs
        self.seq_tech_name = ""
        self.seq_tech_line_edit.textChanged.connect(self.update_seq_tech)
        self.synth_tech_name = ""
        self.synth_tech_line_edit.textChanged.connect(self.update_synth_tech)

        # solqc Qslider initializations
        self.sample_size_label.setVisible(False)
        self.sample_size_slider.setVisible(False)
        self.sample_size_warning_label.setVisible(False)
        self.sample_size_checkbox_state = False
        self.sample_size_slider.setMinimum(10)
        self.sample_size_slider.setMaximum(100)
        self.sample_size_slider.setSingleStep(10)
        self.sample_size_slider.setTickPosition(QSlider.TicksAbove)
        self.slider_val = self.sample_size_slider.value()
        self.sample_size_label.setText(str(self.slider_val) + "%")
        self.sample_size_check_box.clicked.connect(self.choose_sample_size)
        self.sample_size_slider.valueChanged.connect(self.update_slider)

        # set up hints for buttons and labels
        self.design_file_path_label.setToolTip("a csv file containing the design data")
        self.sequencing_file_path_label.setToolTip("a txt file containing the absolute paths to all sequencing (fasta/q) files to be processed")
        self.config_file_path_label.setToolTip("a json file containing the configuration fields")
        self.output_file_path_label.setToolTip("the destination folder for the output file")
        self.seq_tech_label.setToolTip("This field will affect the result file name")
        self.synth_tech_label.setToolTip("This field will affect the result file name")
        self.sample_size_check_box.setToolTip('The percentage selected determines the portion of the sequencing file to be processed')

    def open_design_file_dialog(self):
        self.design_file_path, _ = QFileDialog.getOpenFileName(self, "Select an input file", self.last_dir,
                                                               filter="*.csv")
        self.design_file_path_line_edit.setText(self.design_file_path)
        self.last_dir = self.design_file_path[: self.design_file_path.rfind("/")]

    def open_sequencing_file_dialog(self):
        self.seq_file_path, _ = QFileDialog.getOpenFileName(self, "Select an input file", self.last_dir, filter="*.txt")
        self.sequencing_file_path_line_edit.setText(self.seq_file_path)
        self.last_dir = self.seq_file_path[: self.seq_file_path.rfind("/")]

    def open_config_file_dialog(self):
        self.config_file_path, _ = QFileDialog.getOpenFileName(self, "Select an input file", self.last_dir,
                                                               filter="*.json")
        self.config_file_path_line_edit.setText(self.config_file_path)
        self.last_dir = self.config_file_path[: self.config_file_path.rfind("/")]

    def open_output_folder_dialog(self):
        self.output_folder_path = QFileDialog.getExistingDirectory(self, "Select a folder", self.last_dir)
        self.output_file_path_line_edit.setText(self.output_folder_path)
        self.last_dir = self.output_folder_path

    def search_for_json_file(self):
        json_file_path, _ = QFileDialog.getOpenFileName(self, "Select an input file", self.last_dir, filter="*.json")
        if (json_file_path):
            self.json_path_lineEdit.setText(json_file_path)
            self.last_dir = self.output_folder_path

    def open_and_load_json_file(self):
        self.parseJson(self.json_path_lineEdit.text())

    def parseJson(self, path):
        try:
            with open(path, "r") as read_file:
                tree = json.load(read_file)
        except:
            self.msg_box_with_error("Error opening file. Check path.")
            return
        try:
            self.set_statistics_values(tree['statistics']['substitution'], tree['statistics']['insertion'], tree['statistics']['one_base_deletion'], tree['statistics']['long_deletion'])
            self.set_per_base_substitution(tree['substitution']['a'], tree['substitution']['c'], tree['substitution']['g'], tree['substitution']['t'])
            self.set_per_base_insertion(tree['insert']['a'], tree['insert']['c'], tree['insert']['g'], tree['insert']['t'])
            self.set_per_base_pre_insertion(tree['pre_insert']['a'], tree['pre_insert']['c'], tree['pre_insert']['g'], tree['pre_insert']['t'])
            self.set_per_base_del(tree['one_base_delete']['a'], tree['one_base_delete']['c'], tree['one_base_delete']['g'], tree['one_base_delete']['t'])
            self.set_per_base_long_del(tree['long_delete']['a'], tree['long_delete']['c'], tree['long_delete']['g'], tree['long_delete']['t'])
        except:
            self.msg_box_with_error("Error parsing file. Format might be invalid")

    def msg_box_with_error(self, error_msg):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText(error_msg)
            msg.setWindowTitle("Error!")
            msg.setStandardButtons(QMessageBox.Ok)
            retval = msg.exec_()

    def check_config_file(self) -> bool:
        is_correct = True
        with open(self.config_file_path, "r") as fp:
            data = json.load(fp)
            keys = data.keys()
            keys_list = list(keys)
            keys_list.sort()
            keys_benchmark_list = ["prefix", "suffix", "length", "barcode_start", "barcode_end", "barcode_tolerance"]
            keys_benchmark_list.sort()
            if not keys_list == keys_benchmark_list:
                is_correct = False
            return is_correct

    def check_before_start(self):
        f1 = f2 = f3 = f4 = ""
        error_flag = False
        ######sequence file######
        if not self.design_file_path:
            self.design_file_path_label.setStyleSheet('color: Red')
            f1 = "Design file is missing\n"
            error_flag = True
        elif not self.design_file_path.endswith('.csv'):
            self.design_file_path_label.setStyleSheet('color: Red')
            f1 = "Design file is of wrong type\n"
            error_flag = True
        else:
            self.design_file_path_label.setStyleSheet('color: Black')
        ######sequence file######
        if not self.seq_file_path:
            self.sequencing_file_path_label.setStyleSheet('color: Red')
            f2 = "Sequence file is missing\n"
            error_flag = True
        elif not self.seq_file_path.endswith('.txt'):
            self.sequencing_file_path_label.setStyleSheet('color: Red')
            f2 = "Sequence file is of wrong type\n"
            error_flag = True
        else:
            self.sequencing_file_path_label.setStyleSheet('color: Black')
        ######configuration file######
        if not self.config_file_path:
            self.config_file_path_label.setStyleSheet('color: Red')
            f3 = "Configuration file is missing\n"
            error_flag = True
        elif not self.config_file_path.endswith('.json'):
            self.config_file_path_label.setStyleSheet('color: Red')
            f3 = "Configuration file is of wrong type\n"
            error_flag = True
        elif not self.check_config_file():
            self.config_file_path_label.setStyleSheet('color: Red')
            f3 = "Configuration file is corrupted\n"
            error_flag = True
        else:
            self.config_file_path_label.setStyleSheet('color: Black')
        ######Output Folder######
        if not self.output_folder_path:
            self.output_file_path_label.setStyleSheet('color: Red')
            f4 = "Output Folder is missing"
            error_flag = True
        elif not isdir(self.output_folder_path):
            self.output_file_path_label.setStyleSheet('color: Red')
            f4 = "Output Folder does not exsit"
            error_flag = True
        else:
            self.output_file_path_label.setStyleSheet('color: Black')
        if error_flag == True:
            error_str = f"{f1}{f2}{f3}{f4}"
            self.msg_box_with_error(error_str.strip())
        return error_flag

    def move_temp_files(self):
        # move files from solqc temp folder to selected output folder
        src_dir = os.getcwd() + "/temp"
        dest_dir = self.output_folder_path
        config_name_prefix = ""
        if self.synth_tech_name != "":
            config_name_prefix += self.synth_tech_name + "_synth_"
        if self.seq_tech_name != "":
            config_name_prefix += self.seq_tech_name + "_seq_"
        for file in os.listdir(src_dir):
            if file.find("storalator_config") != -1:
                shutil.move(src_dir + "/" + file, dest_dir + "/" + config_name_prefix + file)
            else:
                os.remove(src_dir + "/" + file)

    def removeTempFiles(self):
        if self.sample_size_checkbox_state == True:
            if os.path.exists(self.output_folder_path + "/temp_data.fastq"):
                os.remove(self.output_folder_path + "/temp_data.fastq")
            if os.path.exists(self.output_folder_path + "/temp_seq.txt"):
                os.remove(self.output_folder_path + "/temp_seq.txt")

    def start_solqc(self):
        # start solqc from command line terminal
        is_error = self.check_before_start()
        if not is_error:
            if self.sample_size_checkbox_state == True and self.slider_val < 100:
                seq_file = self.partial_sample()
                if seq_file == -1:
                    self.removeTempFiles()
                    return
            else:
                seq_file = self.seq_file_path
            cmd = "python solqc/mainStoralator.py -d " + self.design_file_path + " -c " + self.config_file_path + " -r " + seq_file
            try:
                os.system(cmd)
            except:
                self.msg_box_with_error("Error running Solqc")
            self.move_temp_files()
            self.removeTempFiles()

    def partial_sample(self):
        '''returns a path to temporary file contains the percentage of original
        sequencing file that the user chose to process using the slider. the size of the file is the closest
        to percetage specified by user as a complete multiple of 4 to comply with sequence information completeness'''
        f = open(self.seq_file_path, "r")  # path to fastq file
        total_lines = 0
        temp_data_file_path = self.output_folder_path + "/temp_data.fastq"
        for line in f:
            try:
                temp_file = open(line, "r")
            except:
                self.msg_box_with_error("Error opening file: %s" % line)
                return -1
            temp_file = open(line, "r")
            total_lines += len(temp_file.readlines())
            temp_file.close()
        calculated_lines = float(round((total_lines * (self.slider_val/100)) - (total_lines * (self.slider_val/100) % 4)))  # value should come from Qslider value
        try:
            new_file = open(temp_data_file_path,"w+")  # path to new file
            f.seek(0)
            for line in f:
                if calculated_lines == 0:
                    break
                try:
                    temp_file = open(line, "r")
                    while calculated_lines:
                        line = temp_file.readline()
                        print(line)
                        new_file.write(str(line))
                        calculated_lines = calculated_lines - 1
                    temp_file.close()
                except:
                    self.msg_box_with_error("Error opening file: %s" % line)
                    return -1
            new_file.close()
        except:
            self.msg_box_with_error("Error creating partial sample file.")
            return -1
        f.close()
        temp_seq_file_path = self.output_folder_path + "/temp_seq.txt"
        new_file = open(temp_seq_file_path,"w+")  # path to new file
        new_file.write(temp_data_file_path)
        new_file.close()
        return temp_seq_file_path

    def show_solqc_help(self):
        self.solqc_help_window.show()

    def update_seq_tech(self, text):
        self.seq_tech_name = text

    def update_synth_tech(self, text):
        self.synth_tech_name = text

    def choose_sample_size(self):
        if self.sample_size_checkbox_state == False:
            self.sample_size_checkbox_state = True
            self.sample_size_slider.setVisible(True)
            self.sample_size_label.setVisible(True)
            if self.slider_val > 60:
                self.sample_size_warning_label.setVisible(True)
        else:
            self.sample_size_checkbox_state = False
            self.sample_size_slider.setVisible(False)
            self.sample_size_label.setVisible(False)
            self.sample_size_warning_label.setVisible(False)

    def update_slider(self):
        self.slider_val = self.sample_size_slider.value()
        self.sample_size_label.setText(str(self.slider_val) + "%")
        if self.slider_val > 60:
            self.sample_size_warning_label.setVisible(True)
        else:
            self.sample_size_warning_label.setVisible(False)

    def user_defined_amount(self):
        checkbox_status = self.user_defined_copied_checkBox.isChecked()
        self.cont_radioButton.setVisible(checkbox_status)
        self.vector_radioButton.setVisible(checkbox_status)
        if not checkbox_status:
            self.value_lineEdit.setVisible(checkbox_status)
            self.min_amount_spinBox.setVisible(checkbox_status)
            self.max_amount_spinBox.setVisible(checkbox_status)
            self.value_label.setVisible(checkbox_status)
            self.min_label.setVisible(checkbox_status)
            self.max_label.setVisible(checkbox_status)
        elif self.cont_radioButton.isChecked():
            self.continuous_amount()
        elif self.vector_radioButton.isChecked():
            self.vector_amount()

    def continuous_amount(self):
        self.value_lineEdit.setVisible(True)
        self.value_label.setVisible(True)
        self.min_amount_spinBox.setVisible(True)
        self.min_label.setVisible(True)
        self.max_amount_spinBox.setVisible(True)
        self.max_label.setVisible(True)
        self.dist_info['type'] = self.cont_radioButton.text().lower()

    def vector_amount(self):
        self.value_lineEdit.setVisible(True)
        self.value_label.setVisible(True)
        self.min_amount_spinBox.setVisible(False)
        self.min_label.setVisible(False)
        self.max_amount_spinBox.setVisible(False)
        self.max_label.setVisible(False)
        self.dist_info['type'] = self.vector_radioButton.text().lower()

    def set_min_amount(self, value):
        self.dist_info['min'] = int(value)

    def set_max_amount(self, value):
        self.dist_info['max'] = int(value)

    def set_value_amount(self, value):
        self.dist_info['value'] = value

    def set_barcode_start(self, value):
        self.barcode_start = value

    def set_barcode_end(self, value):
        self.barcode_end = value

    def set_clustering_edit_dist(self, value):
        self.max_clustering_edit_dist = value

    def set_minHash_local_steps(self, value):
        self.minHash_local_steps = value
    
    def set_minHash_window_size(self, value):
        self.minHash_window_size = value
    
    def set_minHash_similarity_threshold(self, value):
        self.minHash_similarity_threshold = value

    def set_reconstruction_algo(self):
        self.reconstruction_algo = self.reconstruction_listWidget.currentItem().text()
        # refreshes the view to enable the show button
        self.reconstruction_table.model().dataChanged.emit(QModelIndex(), QModelIndex())

    def set_clustering_algo(self):
        self.clustering_algo = self.clustering_listWidget.currentItem().text()
        self.stackedWidget.show()
        self.clustering_settings_label.setVisible(True)
        clustering_algo_index_dict = {
            'Pseudo Clustering Algorithm': 0,
            'Index Based Algorithm': 1,
            'Hash Based Algorithm': 2
        }
        self.stackedWidget.setCurrentIndex(clustering_algo_index_dict[self.clustering_algo])
        if self.clustering_algo == 'Index Based Algorithm':
            self.set_chosen_clustering_technology()
            self.set_clustering_index()
        # else: # Pseudo, nothing to do

    def set_chosen_clustering_technology(self):
        self.cluster_chosen_technology = self.chosen_tech_dict[self.technology_comboBox.currentText()]
        # refreshes the view to enable the show button
        self.clustering_table.model().dataChanged.emit(QModelIndex(), QModelIndex())
        
    def set_chosen_minHash_clustering_technology(self):
        self.cluster_chosen_technology = self.chosen_tech_dict[self.minHash_technology_comboBox.currentText()]
        # refreshes the view to enable the show button
        self.clustering_table.model().dataChanged.emit(QModelIndex(), QModelIndex())

    def set_clustering_index(self):
        self.clustering_index = int(self.clustering_index_comboBox.currentText())
        # refreshes the view to enable the show button
        self.clustering_table.model().dataChanged.emit(QModelIndex(), QModelIndex())

    def set_substitution(self, value):
        self.general_errors['s'] = value

    def set_insertion(self, value):
        self.general_errors['i'] = value

    def set_one_base_del(self, value):
        self.general_errors['d'] = value

    def set_long_del(self, value):
        self.general_errors['ld'] = value

    def set_A_substitution(self, value):
        self.per_base_errors['A']['s'] = value

    def set_C_substitution(self, value):
        self.per_base_errors['C']['s'] = value

    def set_G_substitution(self, value):
        self.per_base_errors['G']['s'] = value

    def set_T_substitution(self, value):
        self.per_base_errors['T']['s'] = value

    def set_A_insertion(self, value):
        self.per_base_errors['A']['i'] = value

    def set_C_insertion(self, value):
        self.per_base_errors['C']['i'] = value

    def set_G_insertion(self, value):
        self.per_base_errors['G']['i'] = value

    def set_T_insertion(self, value):
        self.per_base_errors['T']['i'] = value

    def set_A_pre_insertion(self, value):
        self.per_base_errors['A']['pi'] = value

    def set_C_pre_insertion(self, value):
        self.per_base_errors['C']['pi'] = value

    def set_G_pre_insertion(self, value):
        self.per_base_errors['G']['pi'] = value

    def set_T_pre_insertion(self, value):
        self.per_base_errors['T']['pi'] = value

    def set_A_one_base_del(self, value):
        self.per_base_errors['A']['d'] = value

    def set_C_one_base_del(self, value):
        self.per_base_errors['C']['d'] = value

    def set_G_one_base_del(self, value):
        self.per_base_errors['G']['d'] = value

    def set_T_one_base_del(self, value):
        self.per_base_errors['T']['d'] = value

    def set_A_long_del(self, value):
        self.per_base_errors['A']['ld'] = value

    def set_C_long_del(self, value):
        self.per_base_errors['C']['ld'] = value

    def set_G_long_del(self, value):
        self.per_base_errors['G']['ld'] = value

    def set_T_long_del(self, value):
        self.per_base_errors['T']['ld'] = value

    def MinIon_chosen(self):
        self.twist_bioscience_radioButton.setEnabled(False)
        self.twist_bioscience_radioButton.setAutoExclusive(False)
        self.twist_bioscience_radioButton.setChecked(False)
        self.twist_bioscience_radioButton.setAutoExclusive(True)
        self.customArray_radioButton.setEnabled(False)
        self.customArray_radioButton.setAutoExclusive(False)
        self.customArray_radioButton.setChecked(False)
        self.customArray_radioButton.setAutoExclusive(True)
        self.IDT_radioButton.setEnabled(True)
        # self.IDT_radioButton.setChecked(True)
        self.default_radioButton.setEnabled(False)
        self.default_radioButton.setAutoExclusive(False)
        self.default_radioButton.setChecked(False)
        self.default_radioButton.setAutoExclusive(True)
        self.error_stats_label.setVisible(True)
        self.error_stats_visibility(True)

    def Stutter_chosen(self):
        self.twist_bioscience_radioButton.setEnabled(False)
        self.twist_bioscience_radioButton.setAutoExclusive(False)
        self.twist_bioscience_radioButton.setChecked(False)
        self.twist_bioscience_radioButton.setAutoExclusive(True)
        self.customArray_radioButton.setEnabled(False)
        self.customArray_radioButton.setAutoExclusive(False)
        self.customArray_radioButton.setChecked(False)
        self.customArray_radioButton.setAutoExclusive(True)
        self.IDT_radioButton.setEnabled(False)
        self.IDT_radioButton.setAutoExclusive(False)
        self.IDT_radioButton.setChecked(False)
        self.IDT_radioButton.setAutoExclusive(True)
        self.default_radioButton.setEnabled(True)

    def error_stats_visibility(self, visibility):
        self.error_stats_label.setVisible(visibility)
        self.substitution_doubleSpinBox.setVisible(visibility)
        self.insertion_doubleSpinBox.setVisible(visibility)
        self.one_base_del_doubleSpinBox.setVisible(visibility)
        self.long_del_doubleSpinBox.setVisible(visibility)
        self.substitution_general_label.setVisible(visibility)
        self.insertion_general_label.setVisible(visibility)
        self.one_base_del_general_label.setVisible(visibility)
        self.long_del_general_label.setVisible(visibility)
        self.line_4.setVisible(visibility)
        self.pre_insertion_label.setVisible(visibility)
        self.A_pre_insertion_doubleSpinBox.setVisible(visibility)
        self.C_pre_insertion_doubleSpinBox.setVisible(visibility)
        self.G_pre_insertion_doubleSpinBox.setVisible(visibility)
        self.T_pre_insertion_doubleSpinBox.setVisible(visibility)
        self.long_del_label.setVisible(visibility)
        self.A_long_del_doubleSpinBox.setVisible(visibility)
        self.C_long_del_doubleSpinBox.setVisible(visibility)
        self.G_long_del_doubleSpinBox.setVisible(visibility)
        self.T_long_del_doubleSpinBox.setVisible(visibility)
        if visibility:
            self.one_base_del_label.setText('1-base Deletion')
            self.Insertion_label.setText('Insertion')
        else:
            self.one_base_del_label.setText('Deletion')
            self.Insertion_label.setText('Stutter')

    def ilumina_NextSeq_chosen(self):
        self.twist_bioscience_radioButton.setEnabled(True)
        # self.twist_bioscience_radioButton.setChecked(True)
        self.customArray_radioButton.setEnabled(False)
        self.customArray_radioButton.setAutoExclusive(False)
        self.customArray_radioButton.setChecked(False)
        self.customArray_radioButton.setAutoExclusive(True)
        self.IDT_radioButton.setEnabled(False)
        self.IDT_radioButton.setAutoExclusive(False)
        self.IDT_radioButton.setChecked(False)
        self.IDT_radioButton.setAutoExclusive(True)
        self.default_radioButton.setEnabled(False)
        self.default_radioButton.setAutoExclusive(False)
        self.default_radioButton.setChecked(False)
        self.default_radioButton.setAutoExclusive(True)
        self.error_stats_label.setVisible(True)
        self.error_stats_visibility(True)

    def ilumina_miSeq_chosen(self):
        self.twist_bioscience_radioButton.setEnabled(True)
        self.customArray_radioButton.setEnabled(True)
        self.IDT_radioButton.setEnabled(False)
        self.IDT_radioButton.setAutoExclusive(False)
        self.IDT_radioButton.setChecked(False)
        self.IDT_radioButton.setAutoExclusive(True)
        self.twist_bioscience_radioButton.setAutoExclusive(False)
        self.twist_bioscience_radioButton.setChecked(False)
        self.twist_bioscience_radioButton.setAutoExclusive(True)
        self.default_radioButton.setEnabled(False)
        self.default_radioButton.setAutoExclusive(False)
        self.default_radioButton.setChecked(False)
        self.default_radioButton.setAutoExclusive(True)
        self.error_stats_label.setVisible(True)
        self.error_stats_visibility(True)

    def twist_bioscience_chosen(self):
        if self.Ilumina_miSeq_radioButton.isChecked() and self.twist_bioscience_radioButton.isChecked():
            self.error_chosen_technology = 'miseq_twist'
            self.set_EZ17_values()
        elif self.Ilumina_NextSeq_radioButton.isChecked() and self.twist_bioscience_radioButton.isChecked():
            self.error_chosen_technology = 'nextseq_twist'
            self.set_O17_values()

    def customArray_chosen(self):
        if self.Ilumina_miSeq_radioButton.isChecked() and self.customArray_radioButton.isChecked():
            self.error_chosen_technology = 'miseq_custom'
            self.set_G15_values()

    def IDT_chosen(self):
        if self.MinION_radioButton.isChecked() and self.IDT_radioButton.isChecked():
            self.error_chosen_technology = 'minion_idt'
            self.set_Y16_values()

    def Default_chosen(self):
        if self.stutter_radioButton.isChecked() and self.default_radioButton.isChecked():
            self.error_chosen_technology = 'stutter_def'
            self.error_stats_visibility(False)
            self.set_stutter_default_values()

    def set_statistics_values(self, a, c, g, t):
        self.substitution_doubleSpinBox.setValue(a)
        self.insertion_doubleSpinBox.setValue(c)
        self.one_base_del_doubleSpinBox.setValue(g)
        self.long_del_doubleSpinBox.setValue(t)

    def set_per_base_substitution(self, a, c, g, t):
        self.A_substitution_doubleSpinBox.setValue(a)
        self.C_substitution_doubleSpinBox.setValue(c)
        self.G_substitution_doubleSpinBox.setValue(g)
        self.T_substitution_doubleSpinBox.setValue(t)

    def set_per_base_insertion(self, a, c, g, t):
        self.A_insertion_doubleSpinBox.setValue(a)
        self.C_insertion_doubleSpinBox.setValue(c)
        self.G_insertion_doubleSpinBox.setValue(g)
        self.T_insertion_doubleSpinBox.setValue(t)

    def set_per_base_pre_insertion(self, a, c, g, t):
        self.A_pre_insertion_doubleSpinBox.setValue(a)
        self.C_pre_insertion_doubleSpinBox.setValue(c)
        self.G_pre_insertion_doubleSpinBox.setValue(g)
        self.T_pre_insertion_doubleSpinBox.setValue(t)

    def set_per_base_del(self, a, c, g, t):
        self.A_one_base_del_doubleSpinBox.setValue(a)
        self.C_one_base_del_doubleSpinBox.setValue(c)
        self.G_one_base_del_doubleSpinBox.setValue(g)
        self.T_one_base_del_doubleSpinBox.setValue(t)

    def set_per_base_long_del(self, a, c, g, t):
        self.A_long_del_doubleSpinBox.setValue(a)
        self.C_long_del_doubleSpinBox.setValue(c)
        self.G_long_del_doubleSpinBox.setValue(g)
        self.T_long_del_doubleSpinBox.setValue(t)

    def set_stutter_default_values(self):
        self.substitution_doubleSpinBox.setValue(0.1)
        self.insertion_doubleSpinBox.setValue(0.1)
        self.one_base_del_doubleSpinBox.setValue(0.1)
        self.long_del_doubleSpinBox.setValue(0.1)

        self.set_per_base_substitution(0.1, 0.1, 0.1, 0.1)
        self.set_per_base_insertion(0.1, 0.1, 0.1, 0.1)
        self.set_per_base_pre_insertion(0.1, 0.1, 0.1, 0.1)
        self.set_per_base_del(0.1, 0.1, 0.1, 0.1)
        self.set_per_base_long_del(0.1, 0.1, 0.1, 0.1)

    def set_minHash_cluster_values(self):
        self.minHash_local_steps_spinBox.setRange(self.minHash_min_local_steps, self.minHash_max_local_steps)
        self.minHash_local_steps_spinBox.setValue(self.minHash_local_steps)
        self.minHash_window_size_spinBox.setRange(self.minHash_min_window_size, self.minHash_max_window_size)
        self.minHash_window_size_spinBox.setValue(self.minHash_window_size)
        self.minHash_similarity_threshold_spinBox.setRange(self.minHash_min_similarity_threshold, self.minHash_max_similarity_threshold)
        self.minHash_similarity_threshold_spinBox.setValue(self.minHash_similarity_threshold)

    def set_EZ17_values(self):
        # general errors
        self.substitution_doubleSpinBox.setValue(1.32e-03)
        self.insertion_doubleSpinBox.setValue(5.81e-04)
        self.one_base_del_doubleSpinBox.setValue(9.58e-04)
        self.long_del_doubleSpinBox.setValue(2.33e-04)

        # per base errors
        self.set_per_base_substitution(0.00135, 0.00135, 0.00126, 0.00132)
        self.set_per_base_insertion(0.00057, 0.00059, 0.00059, 0.00058)
        self.set_per_base_pre_insertion(0.00059, 0.00058, 0.00057, 0.00058)
        self.set_per_base_del(0.00099, 0.00098, 0.00094, 0.00096)
        self.set_per_base_long_del(0.00024, 0.00023, 0.00023, 0.00023)

    def set_O17_values(self):
        # general errors
        self.substitution_doubleSpinBox.setValue(7.09e-03)
        self.insertion_doubleSpinBox.setValue(4.14e-03)
        self.one_base_del_doubleSpinBox.setValue(2.77e-03)
        self.long_del_doubleSpinBox.setValue(4.79e-04)

        # per base errors
        self.set_per_base_substitution(0.00724, 0.00701, 0.00706, 0.00704)
        self.set_per_base_insertion(0.00411, 0.00415, 0.00415, 0.00413)
        self.set_per_base_pre_insertion(0.00429, 0.00415, 0.00403, 0.00408)
        self.set_per_base_del(0.00289, 0.00279, 0.00276, 0.0028)
        self.set_per_base_long_del(0.00048, 0.00048, 0.00047, 0.00049)

    def set_G15_values(self):
        # general errors
        self.substitution_doubleSpinBox.setValue(5.84e-03)
        self.insertion_doubleSpinBox.setValue(8.57e-04)
        self.one_base_del_doubleSpinBox.setValue(5.37e-03)
        self.long_del_doubleSpinBox.setValue(3.48e-04)

        # per base errors
        self.set_per_base_substitution(0.00605, 0.00563, 0.00577, 0.00591)
        self.set_per_base_insertion(0.0009, 0.00083, 0.00085, 0.00084)
        self.set_per_base_pre_insertion(0.00092, 0.00081, 0.00087, 0.00084)
        self.set_per_base_del(0.00543, 0.00513, 0.00539, 0.00559)
        self.set_per_base_long_del(0.00036, 0.00034, 0.00034, 0.00036)

    def set_Y16_values(self):
        # general errors
        self.substitution_doubleSpinBox.setValue(2.2e-02)
        self.insertion_doubleSpinBox.setValue(1.7e-02)
        self.one_base_del_doubleSpinBox.setValue(0.2e-01)
        self.long_del_doubleSpinBox.setValue(0.4e-02)

        # per base errors
        self.set_per_base_substitution(0.119, 0.133, 0.112, 0.119)
        self.set_per_base_insertion(0.331, 0.406, 0.361, 0.367)
        self.set_per_base_pre_insertion(0.332, 0.408, 0.341, 0.382)
        self.set_per_base_del(0.044, 0.048, 0.040, 0.041)
        self.set_per_base_long_del(0.019, 0.021, 0.017, 0.018)

    def openFileDialog(self):
        self.inputDNAPath, _ = QFileDialog.getOpenFileName(self, "Select an input file", './', filter="*.txt")
        self.file_path_lineEdit.setText(self.inputDNAPath)

    def set_working_dirs(self):
        self.input_dir = "input"
        self.output_dir = "output"

        self.error_sim_output_dir = self.output_dir + '/error_sim'
        self.clustering_output_dir = self.output_dir + '/clustering'
        self.clustering_pseudo_dir = self.clustering_output_dir + '/pseudo'
        self.clustering_index_dir = self.clustering_output_dir + '/index'
        self.clustering_hash_dir = self.clustering_output_dir + '/hash'
        self.clustering_result_dir = self.clustering_output_dir + '/res'

        self.recon_output_dir = self.output_dir + '/reconstruction'
        self.recon_algo_dir = 'reconstruction_algs'
        if platform.system() == "Windows":
            self.recon_algo_file_suffix = '.exe'
        elif platform.system() == "Linux":
            self.recon_algo_file_suffix = 'Linux'
        elif platform.system() == "Darwin":
            self.recon_algo_file_suffix = 'Mac'
        self.recon_hist_dir = self.recon_output_dir + "/histograms"

    def create_base_dir(self):
        os.mkdir(self.output_dir)

    def set_json(self):
        self.recon_json = self.recon_output_dir + '/recon_results.json'
        self.clustering_res_json = self.clustering_result_dir + '/results.json'

    def error_sim_output_files_path(self, timestamp):
        suffix = '_' + self.error_chosen_technology + '_' + str(timestamp) + '.txt'
        error_path, shuffled_path = self.error_sim_output_dir + '/error' + suffix, self.error_sim_output_dir + '/error_shuffled' + suffix
        return error_path, shuffled_path

    def clustering_output_file_path(self, timestamp):
        if self.clustering_algo == 'Pseudo Clustering Algorithm':
            dir = '/pseudo'
            index = ''
        elif self.clustering_algo == 'Hash Based Algorithm':
            dir = '/hash'
            index = ''
        else:
            dir = '/index'
            index = str(self.clustering_index)

        suffix = '_' + self.cluster_chosen_technology
        if index != '':
            suffix += '_' + index
        suffix += '_' + str(timestamp) + '.txt'
        return self.clustering_output_dir + dir + '/clustered' + suffix

    def recon_output_file_path(self, timestamp):
        suffix = '-' + str(self.recon_algo) + '_' + str(timestamp) + '.txt'
        return self.recon_output_dir + '/recon' + suffix

    def runErrorSimulator(self):
        self.inputDNAPath = self.file_path_lineEdit.text()
        while True:
            if not os.path.isfile(self.inputDNAPath):
                print('The chosen input file doesn\'t exist')
                self.msg_box_with_title("The input file you chosen doesn't exist", "Error!")
                self.file_path_lineEdit.clear()
                break
            else:
                if self.user_defined_copied_checkBox.isChecked():
                    sent_distance_info = self.dist_info
                    if sent_distance_info['type'] == 'continuous' and sent_distance_info['min'] >= sent_distance_info['max']:
                        self.msg_box_with_title('Min should be < Max in user defined amount of copies', 'Error!')
                        break
                    if sent_distance_info['value'] == '':
                        self.msg_box_with_title('Please enter a value in user defined amount on copies', 'Error!')
                        break
                else:
                    sent_distance_info = None
                # construct paths for output files based on chosen technology and current time:
                self.error_output, self.shuffled_output = self.error_sim_output_files_path(time_stamp())

                error_worker = SimulateErrorsWorker(self.general_errors, self.per_base_errors, self.inputDNAPath,
                                                    self.default_radioButton.isChecked(), sent_distance_info,
                                                    self.error_output, self.shuffled_output)
                # all currently active workers are stored at one place for management:
                self.error_workers.append(error_worker)
                # connect relevant signals:
                error_worker.work_done.connect(lambda: self.error_sim_finished(error_worker))
                error_worker.terminated.connect(lambda: self.error_sim_terminated(error_worker))

                self.error_sim_table_add_file(self.error_output, error_worker, True)
                error_worker.start()
                break

    def error_sim_finished(self, error_worker):
        # update size and timestamp in a relevant row of the table:
        filepath = error_worker.input_path
        self.error_sim_table_update(filepath)
        # new output is a legit input for clustering and reconstruction:
        self.clustering_table_add_file(filepath)
        self.recon_table_add_file(filepath)
        # release worker:
        self.error_workers.remove(error_worker)
        # update progress to 100:
        self.error_sim_table.model().update_cell('100', filepath, self.error_sim_table_cols['progress'])

    def error_sim_terminated(self, error_worker):
        # delete output files, including temp files:
        TableModel.delete_file(self.tables, error_worker.input_path)
        TableModel.delete_file(None, error_worker.shuffled_output_path)
        TableModel.delete_file(None, error_worker.temp_file)
        # release worker:
        self.error_workers.remove(error_worker)

    def evyat_files_path(self):
        if platform.system() == "Linux":
            evyat_path = 'files/' + self.chosen_technology + '/' + 'evyat.txt'
            shuffled_path = 'files/' + self.chosen_technology + '/' + 'errors_shuffled.txt'
            if not os.access('files', os.F_OK):
                os.mkdir('files')
            if not os.access('files/' + self.chosen_technology, os.F_OK):
                os.mkdir('files/' + self.chosen_technology)

        elif platform.system() == "Windows":
            # inputDNAPath = 'files/' + str(index) + '_allclustersofindex.txt'
            evyat_path = 'files/' + self.chosen_technology + '/' + 'evyat.txt'
            shuffled_path = 'files/' + self.chosen_technology + '/' + 'errors_shuffled.txt'
            if not os.access('files', os.F_OK):
                os.mkdir('files')
            if not os.access('files/' + self.chosen_technology, os.F_OK):
                os.mkdir('files/' + self.chosen_technology)
        elif platform.system() == "Darwin":
            # inputDNAPath = 'files/' + str(index) + '_allclustersofindex.txt'
            evyat_path = 'files/' + self.chosen_technology + '/' + 'evyat.txt'
            shuffled_path = 'files/' + self.chosen_technology + '/' + 'errors_shuffled.txt'
            if not os.access('files', os.F_OK):
                os.mkdir('files')
            if not os.access('files/' + self.chosen_technology, os.F_OK):
                os.mkdir('files/' + self.chosen_technology)
    
        return evyat_path, shuffled_path

    def runClustering(self):
        if self.clustering_orig_input == '':
            self.msg_box_with_title("Please, choose an input for Clustering", "Error!")
            return

        # clustering input is set by selecting the row.
        self.clustering_output = self.clustering_output_file_path(time_stamp())
        if self.clustering_algo == 'Pseudo Clustering Algorithm':
            pseudo_cluster(self.clustering_orig_input, self.clustering_output,
                           int(self.barcode_start), int(self.barcode_end), int(self.max_clustering_edit_dist))
            self.recon_table_add_file(self.clustering_output)
            self.msg_box_with_title("Pseudo Clustering successfully completed \n"
                                    "Output file added to Reconstruction tab", "Message!")
        elif self.clustering_algo == 'Index Based Algorithm':

            cluster_worker = ClusteringWorker(self.cluster_chosen_technology, self.clustering_index,
                                              self.clustering_result_dir,
                                              self.clustering_shuffled_input, self.clustering_orig_input,
                                              self.clustering_output)

            # success will hold True iff there is no currently active clustering worker with the same parameters:
            success = self.clustering_table_add_progress(self.clustering_orig_input, self.cluster_chosen_technology,
                                                         str(self.clustering_index), cluster_worker)
            if success:
                self.cluster_workers.append(cluster_worker)
                cluster_worker.work_done.connect(lambda: self.cluster_worker_finished(cluster_worker))
                cluster_worker.terminated.connect(lambda: self.cluster_worker_terminated(cluster_worker))
                cluster_worker.start()
            else:
                self.msg_box_with_title("Already in progress", "Error!")
        elif self.clustering_algo == 'Hash Based Algorithm':
            
            # self.chosen_index_save = self.clustering_index
            self.clustering_index = ""
            #self.evyat_path, self.shuffled_path = self.evyat_files_path()
            cluster_worker = HashBasedClusteringWorker(self.cluster_chosen_technology,self.clustering_index, int(self.minHash_local_steps), int(self.minHash_window_size), int(self.minHash_similarity_threshold), self.clustering_algo, self.clustering_result_dir, self.clustering_orig_input,self.clustering_output)
            # self.cluster_resutls_file = self.cluster_worker.get_results_file()
            # self.progressBar.setValue(0)
            # self.cluster_worker.start()
            # self.label_progress.setText('Clustering in progress, please wait!')
            # self.cluster_worker.finished.connect(self.evt_cluster_worker_finished)
            # self.cluster_worker.update_progress.connect(self.evt_update_progress)
            # self.progressBar.setVisible(True)
            
            success = self.clustering_table_add_progress(self.clustering_orig_input, self.cluster_chosen_technology,
                                                         str(self.clustering_index), cluster_worker)
            if success:
                self.cluster_workers.append(cluster_worker)
                cluster_worker.work_done.connect(lambda: self.cluster_worker_finished(cluster_worker))
                cluster_worker.terminated.connect(lambda: self.cluster_worker_terminated(cluster_worker))
                cluster_worker.start()
            else:
                self.msg_box_with_title("Already in progress", "Error!")

    def cluster_worker_finished(self, cluster_worker):
        # release worker:
        self.cluster_workers.remove(cluster_worker)

        filepath = cluster_worker.input_path
        search_info = dict(data=[cluster_worker.cluster_tech, cluster_worker.cluster_index],
                           cols=[self.cluster_progress_cols['tech'],
                                 self.cluster_progress_cols['index']])
      
        # disable stop button:
        self.clustering_progress_table.model().update_cell('disabled', filepath, self.cluster_progress_cols['stop'],
                                                           in_root=True, secondary_info=search_info, role=Qt.UserRole)
        # update progress to 100:
        self.clustering_progress_table.model().update_cell('100', filepath, self.cluster_progress_cols['progress'],
                                                           in_root=True, secondary_info=search_info)
        # refreshes the view to enable the show button
        self.clustering_table.model().dataChanged.emit(QModelIndex(), QModelIndex())
        # add clustered result to reconstruction table:
        # TODO: since index based clustering doesn't produce an output file, as for now, this feature should be added in the future.
        #self.recon_table_add_file(self.clustering_output)

    def cluster_worker_terminated(self, cluster_worker):
        # release worker:
        self.cluster_workers.remove(cluster_worker)
        # remove row from progress table:
        path = cluster_worker.input_path
        tech = cluster_worker.cluster_tech
        index = str(cluster_worker.cluster_index)
        self.clustering_progress_table.model().remove_row(path, True, dict(data=[tech, index],
                                                                           cols=[self.cluster_progress_cols['tech'],
                                                                                 self.cluster_progress_cols['index']]))
        # clustering intermediate computation files:
        to_delete = [Clustering.outputdict_path(path, tech, index), Clustering.outputstringsdict_path(path, tech, index),
                     Clustering.orig_dict_path(path, tech, index), Clustering.strands_in_wrong_cluster(path, tech, index),
                     Clustering.false_positive(path, tech, index), Clustering.edit_dist_avg_f(path, tech, index)]

        for file in to_delete:
            TableModel.delete_file(None, file)

    def run_reconstruction_algo(self):
        if self.recon_input == '':
            self.msg_box_with_title("Please, choose an input for Reconstruction", "Error!")
            return

        if self.reconstruction_algo not in self.recon_algo_dict:
            self.msg_box_with_title('Please choose a reconstruction algorithm', "Error!")
            return
        self.recon_algo = self.recon_algo_dict[self.reconstruction_algo]
        recon_output = self.recon_output_file_path(time_stamp())

        if platform.system() == "Windows":
            process_path = self.recon_algo_dir + "\\" + self.recon_algo + self.recon_algo_file_suffix
        else:
            process_path = "./" + self.recon_algo_dir + "/" + self.recon_algo + self.recon_algo_file_suffix

        arguments = [self.recon_input, self.recon_output_dir]
        assert (os.path.exists(self.recon_input))
        assert (os.path.exists(process_path))
        with open(self.recon_input, 'r') as input_file:
            assert (input_file.read() != '')
        input_file.close()


        process = ReconProcess(self.recon_input, self.recon_algo, self)
        # not allowing multiple simultaneous runs that are the same:
        success = self.recon_table_add_progress(self.recon_input, self.recon_algo, process, recon_output)
        if success:
            self.recon_processes.append(process)
            process.start(process_path, arguments)

            process.finished.connect(lambda _, exit_status:
                                     self.reconstruction_finished(process, self.recon_input, recon_output, exit_status))
        else:
            self.msg_box_with_title("Already in progress", "Error!")

    def recon_dataReady(self, model, name, algo, process, recon_output_file_path):
        raw_line = str(process.readAll(), 'utf-8')
        line = re.split('\r\n', raw_line.strip())
        for i in line:
            if i.isnumeric():
                model.update_cell(i,
                                  name,
                                  self.recon_progress_cols['progress'],
                                  True,
                                  dict(data=[algo],
                                       cols=[self.recon_progress_cols['algo']]))
        with open(recon_output_file_path, 'a') as output_file:
            output_file.write(raw_line + "\n")

    def reconstruction_finished(self, process, recon_input_file_path, recon_output_file_path, exit_status):
        filepath = recon_input_file_path
        algo = process.algo
        search_info = dict(data=[algo], cols=[self.recon_progress_cols['algo']])

        if exit_status == QProcess.NormalExit:
            with open(recon_output_file_path, 'r') as recon_output_file:
                recon_raw_output = recon_output_file.read()
            assert (recon_raw_output != '')
            json_recon_result = parse_recon_output(self.recon_hist_dir, filepath, recon_output_file_path,
                                                   algo=algo, raw_output=recon_raw_output)

            keys = dict(input_file=filepath, algo=algo)
            update_json_file(self.recon_json, json_recon_result.__dict__, keys)

            # set progress to 100
            self.recon_progress_table.model().update_cell('100', filepath, self.recon_progress_cols['progress'],
                                                          in_root=True, secondary_info=search_info)
            # disable stop button:
            self.recon_progress_table.model().update_cell('disabled', filepath, self.recon_progress_cols['stop'],
                                                          in_root=True, secondary_info=search_info, role=Qt.UserRole)

        else:  # PROCESS CRASHED or terminated
            # remove row from progress table:
            self.recon_progress_table.model().remove_row(filepath, True, search_info)
        # we dont need the output file, only the results in json.
        delete(recon_output_file_path)
        self.recon_processes.remove(process)

    def init_tables(self):
        # init main tables that represent the datasets:
        for table, cols in zip(self.tables, self.table_cols):
            self.init_table(table=table,
                            cols=cols,
                            groups=self.chosen_tech_dict,
                            expandable=True)
        # init helper tables that show the progress of clustering/reconstruction workers:
        for table, cols in zip(self.progress_tables, self.progress_table_cols):
            self.init_table(table=table,
                            cols=cols,
                            groups=None,
                            expandable=False)

    def init_table(self, table, cols, groups=None, expandable=False):
        """
        Here the main work for setting up settings for a table is done like: setting the model, headers, delegates, proxies.
        """
        model = TableModel(cols.keys(), table)
        table.setExpandsOnDoubleClick(False)
        table.setRootIsDecorated(False)
        table.header().setDefaultAlignment(Qt.AlignCenter)

        # set proxies for filtering:
        if table in [self.clustering_table, self.reconstruction_table]:
            # will set the base model as well:
            proxy = TableModel.set_proxy(table, model, cols)
            table.selectionModel().currentRowChanged.connect(lambda cur, prev: self.select_input(table, cur))
            if table == self.clustering_table:
                self.clustering_filter.textChanged.connect(lambda text: FilterModel.filter(proxy, text))
            else:  # reconstruction table:
                self.recon_filter.textChanged.connect(lambda text: FilterModel.filter(proxy, text))
        else:  # other tables don't support filtering:
            table.setModel(model)

        if expandable:
            table.clicked.connect(lambda index: model.on_clicked(index, table, expandable))
            delegate_icon = IconRenameDelegate(self, cols['name'], table)
            table.setItemDelegateForColumn(0, delegate_icon)
        if 'action' in cols:
            delegate_delete = DeleteButtonDelegate(ui=self, col=cols['action'], parent=table)
            table.setItemDelegateForColumn(cols['action'], delegate_delete)
        if 'result' in cols:
            if table == self.clustering_table:
                delegate_show = ClusterShowButtonDelegate(ui=self, col=cols['result'], parent=table)
            else:  # reconstruction table
                delegate_show = ReconShowButtonDelegate(ui=self, col=cols['result'], parent=table)
            table.setItemDelegateForColumn(cols['result'], delegate_show)
        if 'progress' in cols:
            delegate_progress = ProgressBarDelegate(table, expandable)
            table.setItemDelegateForColumn(cols['progress'], delegate_progress)
        if 'stop' in cols:
            delegate_stop = StopButtonDelegate(ui=self, col=cols['stop'], parent=table)
            table.setItemDelegateForColumn(cols['stop'], delegate_stop)

        table.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)

        # set expandable rows in groups:
        if groups is not None:
            for group in groups:
                model.add_group(group)

    def fill_tables(self):
        # tables and dirs they load files from:
        for table, dirs in zip(self.tables, [[self.error_sim_output_dir],
                                             [self.error_sim_output_dir],
                                             [self.error_sim_output_dir, self.clustering_pseudo_dir,
                                              self.clustering_index_dir]]):
            self.fill_table(table, dirs)

    def fill_table(self, table, dirs):
        for dir in dirs:
            for filename in files(dir):
                if 'shuffled' not in filename:
                    filepath = join_path(dir, filename)
                    if table == self.error_sim_table:
                        self.error_sim_table_add_file(filepath)
                    if table == self.clustering_table:
                        self.clustering_table_add_file(filepath)
                    if table == self.reconstruction_table:
                        self.recon_table_add_file(filepath)

    def error_sim_table_add_file(self, filepath, error_worker=None, toExpand=False):
        dummy = '              '
        model = self.error_sim_table.model()
        row_elements = {self.error_sim_table_cols['name']: CellData(get_show_name(filepath), filepath).__dict__,
                        self.error_sim_table_cols['action']: CellData('       delete       ').__dict__}

        if error_worker is not None:  # creation of the dataset is in progress, size/time is not available:
            row_elements[self.error_sim_table_cols['progress']] = CellData('0').__dict__
            error_worker.update_progress.connect(lambda x: model.update_cell(x,
                                                                             filepath,
                                                                             self.error_sim_table_cols['progress']))
            row_elements[self.error_sim_table_cols['size']] = CellData(dummy).__dict__
            row_elements[self.error_sim_table_cols['time']] = CellData(dummy).__dict__

        else:  # dataset exist therefore size/time available.
            row_elements[self.error_sim_table_cols['size']] = CellData(sizeof_fmt(os.path.getsize(filepath))).__dict__
            row_elements[self.error_sim_table_cols['progress']] = CellData('100').__dict__
            creation_stamp = round(os.path.getmtime(filepath))
            row_elements[self.error_sim_table_cols['time']] = CellData(
                str(datetime.fromtimestamp(creation_stamp))).__dict__

        model.append_row_to_group(row_elements, None, toExpand)

    def error_sim_table_update(self, filepath):
        # update size and time of creation for the relevant dataset that was produced:
        model = self.error_sim_table.model()
        size = sizeof_fmt(os.path.getsize(filepath))
        creation_stamp = round(os.path.getmtime(filepath))
        model.update_cell(str(datetime.fromtimestamp(creation_stamp)), filepath, self.error_sim_table_cols['time'])
        model.update_cell(size, filepath, self.error_sim_table_cols['size'])

    def clustering_table_add_file(self, filepath):
        model = TableModel.get_source_model(self.clustering_table)
        row_elements = {self.cluster_table_cols['name']: CellData(get_show_name(filepath), filepath).__dict__,
                        self.cluster_table_cols['action']: CellData('       delete       ').__dict__,
                        self.cluster_table_cols['size']: CellData(sizeof_fmt(os.path.getsize(filepath))).__dict__,
                        self.cluster_table_cols['result']: CellData('       show       ').__dict__,
                        self.cluster_table_cols['']: CellData('').__dict__}

        creation_stamp = round(os.path.getmtime(filepath))
        row_elements[self.cluster_table_cols['time']] = CellData(str(datetime.fromtimestamp(creation_stamp))).__dict__

        model.append_row_to_group(row_elements)

    def recon_table_add_file(self, filepath):
        model = TableModel.get_source_model(self.reconstruction_table)
        row_elements = {self.recon_table_cols['name']: CellData(get_show_name(filepath), filepath).__dict__,
                        self.recon_table_cols['action']: CellData('       delete       ').__dict__,
                        self.recon_table_cols['size']: CellData(sizeof_fmt(os.path.getsize(filepath))).__dict__,
                        self.recon_table_cols['result']: CellData('       show       ').__dict__,
                        self.recon_table_cols['']: CellData('').__dict__}

        creation_stamp = round(os.path.getmtime(filepath))
        row_elements[self.recon_table_cols['time']] = CellData(str(datetime.fromtimestamp(creation_stamp))).__dict__

        model.append_row_to_group(row_elements)

    def clustering_table_add_progress(self, name, tech, index, cluster_worker):
        model = self.clustering_progress_table.model()
        row_elements = {self.cluster_progress_cols['name']: CellData(get_base_name(name), name).__dict__,
                        self.cluster_progress_cols['tech']: CellData(tech).__dict__,
                        self.cluster_progress_cols['index']: CellData(index).__dict__,
                        self.cluster_progress_cols['stop']: CellData('       cancel       ', 'enabled').__dict__,
                        self.cluster_progress_cols['progress']: CellData('0').__dict__}

        search_info = dict(data=[tech, index],
                           cols=[self.cluster_progress_cols['tech'],
                                 self.cluster_progress_cols['index']])

        cluster_worker.update_progress.connect(lambda x: model.update_cell(x,
                                                                           name,
                                                                           self.cluster_progress_cols['progress'],
                                                                           True,
                                                                           search_info))
        return model.append_row_to_group(row_elements, search_info)

    def recon_table_add_progress(self, name, algo, process, recon_output_file):
        model = self.recon_progress_table.model()
        row_elements = {self.recon_progress_cols['name']: CellData(get_base_name(name), name).__dict__,
                        self.recon_progress_cols['algo']: CellData(algo).__dict__,
                        self.recon_progress_cols['stop']: CellData('       cancel       ', 'enabled').__dict__,
                        self.recon_progress_cols['progress']: CellData('0').__dict__}

        search_info = dict(data=[algo],
                           cols=[self.recon_progress_cols['algo']])

        process.readyRead.connect(lambda: self.recon_dataReady(model, name, algo, process, recon_output_file))
        return model.append_row_to_group(row_elements, search_info)

    def select_input(self, table, index):
        if index.parent().isValid():
            filepath = index.siblingAtColumn(0).data(Qt.UserRole)
            if table == self.clustering_table:
                self.clustering_orig_input = filepath
                self.clustering_shuffled_input = get_shuffled_path_from_orig(filepath)
            if table == self.reconstruction_table:
                self.recon_input = filepath

    def msg_box_with_title(self, error_msg, title):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(error_msg)
        msg.setWindowTitle(title)
        msg.setStandardButtons(QMessageBox.Ok)
        retval = msg.exec_()

# class HashBasedClusteringWorker(QThread):
#     update_progress = pyqtSignal(int)

#     def __init__(self, cluster_tech, number_of_steps, windows_size, similarity_threshold, chosen_clustering_algo, clustering_result_dir, clustering_input, clustering_output):
#         super(HashBasedClusteringWorker, self).__init__()
#         self.cluster_tech = cluster_tech
#         self.chosen_clustering_algo = chosen_clustering_algo
#         self.number_of_steps = number_of_steps
#         self.windows_size = windows_size
#         self.similarity_threshold = similarity_threshold
#         self.results_f = clustering_result_dir + '_final_results.txt'
#         self.clustering_input = clustering_input
#         self.clustering_output = clustering_output
#         # self.directory = 'cluster_output/' + self.cluster_tech + '/' + self.chosen_clustering_algo
#         # if not os.path.exists(self.directory):
#         #     os.makedirs(self.directory)	

#     def report_func(self, total_lines, curr_line):
#         percent = int(curr_line * 100 // total_lines)
#         self.update_progress.emit(percent)

#     def get_results_file(self):
#         return self.results_f

#     def run(self):
#         clustering_alg = HashBasedCluster(self.cluster_tech)
#         start_time = time.time()
#         clustering_alg.hash_based_cluster(self.number_of_steps,  self.windows_size,
#                                    self.similarity_threshold, self.report_func,self.clustering_input,self.clustering_output)

#         total_clusters_count = clustering_alg.get_number_of_clusters()
#         with open(self.results_f, 'w') as results_file:
#             print("Run time: %s seconds" % round((time.time() - start_time),2), file=results_file)
#             print('Number of clusters generated: ' + str(total_clusters_count), file=results_file)

class solqcHelpWindow(QMainWindow):
    def __init__(self) -> None:
        super(solqcHelpWindow, self).__init__()
        self.ui = solqc_help_window.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.exit_button.clicked.connect(self.close)

if __name__ == '__main__':
    # Create the Qt Application
    app = QApplication(sys.argv)
    # Create and show the form
    dnaSimulator = dnaSimulator()
    dnaSimulator.show()
    # Run the main Qt loop
    sys.exit(app.exec_())
