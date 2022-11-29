import os
import re
import shlex
import subprocess
import time
from functools import partial
import ast

# from PIL.ImageQt import ImageQt
from PyQt5 import QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import *
import sys
import platform

import dnaSimulator_ui2
# from SpinBoxCustom import SpinBoxCustom
from simulator import *
from index_clustering import *


class dnaSimulator(QMainWindow, dnaSimulator_ui2.Ui_dnaSimulator):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)

        # set the title
        self.setWindowTitle('DNA Simulator')

        self.inputDNAPath = ""
        # self.reconstruction_input_path =

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

        self.barcode_start = 0
        self.barcode_end = 0
        self.max_clustering_edit_dist = 0

        self.reconstruction_algo = ''
        self.clustering_algo = ''
        self.chosen_technology = 'miseq_twist'
        self.clustering_index = 4

        self.evyat_path = ''
        self.shuffled_path = ''

        self.progressBar.setVisible(False)

        # connect push buttons to an event
        self.browse_PushButton.clicked.connect(self.openFileDialog)
        # self.set_current_values_PushButton.clicked.connect(self.setErrorValues)
        self.run_error_simulator_PushButton.clicked.connect(self.runErrorSimulator)
        self.reconstruction_run_pushButton.clicked.connect(self.run_reconstruction_algo)
        self.run_clustering_pushButton.clicked.connect(self.runClustering)

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

        self.reconstruction_listWidget.addItem('Hybrid Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Divider BMA Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('BMA Look Ahead Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Iterative Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Mock Reconstruction Algorithm')
        self.reconstruction_listWidget.addItem('Not a real algo')

        self.reconstruction_listWidget.currentItemChanged.connect(self.set_reconstruction_algo)

        self.clustering_listWidget.addItem('Pseudo Clustering Algorithm')
        self.clustering_listWidget.addItem('Index Based Algorithm')

        self.stackedWidget.addWidget(self.pseudoclustering_settings_verticalWidget)
        self.stackedWidget.addWidget(self.realclustering_settings_verticalWidget)
        self.stackedWidget.hide()
        self.clustering_settings_label.setVisible(False)

        self.technology_comboBox.addItem('Twist Bioscience + Ilumina miSeq')
        self.technology_comboBox.addItem('CustomArray + Ilumina miSeq')
        self.technology_comboBox.addItem('Twist Bioscience + Ilumina NextSeq')
        self.technology_comboBox.addItem('Integrated DNA Technology (IDT) + MinION')
        self.technology_comboBox.currentTextChanged.connect(self.set_chosen_clustering_technology)

        self.reconstrcution_technology_comboBox.addItem('Twist Bioscience + Ilumina miSeq')
        self.reconstrcution_technology_comboBox.addItem('CustomArray + Ilumina miSeq')
        self.reconstrcution_technology_comboBox.addItem('Twist Bioscience + Ilumina NextSeq')
        self.reconstrcution_technology_comboBox.addItem('Integrated DNA Technology (IDT) + MinION')
        self.reconstrcution_technology_comboBox.currentTextChanged.connect(self.set_chosen_reconstruction_technology)

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

    def set_reconstruction_algo(self):
        self.reconstruction_algo = self.reconstruction_listWidget.currentItem().text()
        print(self.reconstruction_algo)

    def set_clustering_algo(self):
        self.clustering_algo = self.clustering_listWidget.currentItem().text()
        print(self.clustering_algo)
        self.stackedWidget.show()
        self.clustering_settings_label.setVisible(True)
        if self.clustering_algo == 'Pseudo Clustering Algorithm':
            self.stackedWidget.setCurrentIndex(0)
        elif self.clustering_algo == 'Index Based Algorithm':
            self.stackedWidget.setCurrentIndex(1)

    def set_chosen_clustering_technology(self):
        clustering_tech = self.technology_comboBox.currentText()
        if clustering_tech == 'Twist Bioscience + Ilumina miSeq':
            self.chosen_technology = 'miseq_twist'
        elif clustering_tech == 'CustomArray + Ilumina miSeq':
            self.chosen_technology = 'miseq_custom'
        elif clustering_tech == 'Twist Bioscience + Ilumina NextSeq':
            self.chosen_technology = 'nextseq_twist'
        elif clustering_tech == 'Integrated DNA Technology (IDT) + MinION':
            self.chosen_technology = 'minion_idt'

    def set_chosen_reconstruction_technology(self):
        reconstruction_tech = self.reconstrcution_technology_comboBox.currentText()
        if reconstruction_tech == 'Twist Bioscience + Ilumina miSeq':
            self.chosen_technology = 'miseq_twist'
        elif reconstruction_tech == 'CustomArray + Ilumina miSeq':
            self.chosen_technology = 'miseq_custom'
        elif reconstruction_tech == 'Twist Bioscience + Ilumina NextSeq':
            self.chosen_technology = 'nextseq_twist'
        elif reconstruction_tech == 'Integrated DNA Technology (IDT) + MinION':
            self.chosen_technology = 'minion_idt'

    def set_clustering_index(self):
        self.clustering_index = int(self.clustering_index_comboBox.currentText())

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
            self.chosen_technology = 'miseq_twist'
            self.set_EZ17_values()
        elif self.Ilumina_NextSeq_radioButton.isChecked() and self.twist_bioscience_radioButton.isChecked():
            self.chosen_technology = 'nextseq_twist'
            self.set_O17_values()

    def customArray_chosen(self):
        if self.Ilumina_miSeq_radioButton.isChecked() and self.customArray_radioButton.isChecked():
            self.chosen_technology = 'miseq_custom'
            self.set_G15_values()

    def IDT_chosen(self):
        if self.MinION_radioButton.isChecked() and self.IDT_radioButton.isChecked():
            self.chosen_technology = 'minion_idt'
            self.set_Y16_values()

    def Default_chosen(self):
        if self.stutter_radioButton.isChecked() and self.default_radioButton.isChecked():
            self.error_stats_visibility(False)
            self.set_stutter_default_values()

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
        self.substitution_doubleSpinBox.setValue(1.21e-01)
        self.insertion_doubleSpinBox.setValue(3.67e-01)
        self.one_base_del_doubleSpinBox.setValue(4.33e-02)
        self.long_del_doubleSpinBox.setValue(1.87e-02)

        # per base errors
        self.set_per_base_substitution(0.119, 0.133, 0.112, 0.119)
        self.set_per_base_insertion(0.331, 0.406, 0.361, 0.367)
        self.set_per_base_pre_insertion(0.332, 0.408, 0.341, 0.382)
        self.set_per_base_del(0.044, 0.048, 0.040, 0.041)
        self.set_per_base_long_del(0.019, 0.021, 0.017, 0.018)

    def openFileDialog(self):
        self.inputDNAPath, _ = QFileDialog.getOpenFileName(self, "Select an input file", './', filter="*.txt")
        self.file_path_lineEdit.setText(self.inputDNAPath)

    # def show_error_dialog(self, error_msg):
    #     msg = QMessageBox()
    #     msg.setIcon(QMessageBox.Information)
    #
    #     msg.setText(error_msg)
    #     # msg.setInformativeText("This is additional information")
    #     msg.setWindowTitle("Error!")
    #     # msg.setDetailedText("The details are as follows:")
    #     msg.setStandardButtons(QMessageBox.Ok)
    #     # msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
    #     # msg.buttonClicked.connect(msgbtn)

    # retval = msg.exec_()

    def evyat_files_path(self):
        if platform.system() == "Linux":
            # inputDNAPath = '/home_nfs/sgamer8/DNAindex' + str(index) + '/files/' + str(index) + '_allclustersofindex.txt'
            evyat_path = '/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology + '/' + str(
                self.clustering_index) + '_evyat.txt'
            shuffled_path = '/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology + '/' + str(
                self.clustering_index) + '_errors_shuffled.txt'
            if not os.access('/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology, os.F_OK):
                os.mkdir('/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology)
        elif platform.system() == "Windows":
            # inputDNAPath = 'files/' + str(index) + '_allclustersofindex.txt'
            evyat_path = 'files/' + self.chosen_technology + '/' + 'evyat.txt'
            shuffled_path = 'files/' + self.chosen_technology + '/' + 'errors_shuffled.txt'
            if not os.access('files', os.F_OK):
                os.mkdir('files')
            if not os.access('files/' + self.chosen_technology, os.F_OK):
                os.mkdir('files/' + self.chosen_technology)

        return evyat_path, shuffled_path

    def runClustering(self):
        if self.clustering_algo == 'Pseudo Clustering Algorithm':
            if not os.path.isfile('output/evyat.txt'):
                print('The evyat.txt doesn\'t exist')
                self.msg_box_with_error("evyat.txt input file for clustering doesn't exist")
            else:
                pseudo_cluster(int(self.barcode_start), int(self.barcode_end), int(self.max_clustering_edit_dist))
        elif self.clustering_algo == 'Index Based Algorithm':
            print(self.chosen_technology)
            self.evyat_path, self.shuffled_path = self.evyat_files_path()
            self.cluster_worker = ClusteringWorker(self.chosen_technology, self.clustering_index)
            self.progressBar.setValue(0)
            self.cluster_worker.start()
            self.label_progress.setText('Clustering in progress, please wait!')
            self.clustering_results_textEdit.clear()
            self.cluster_worker.finished.connect(self.evt_cluster_worker_finished)
            self.cluster_worker.update_progress.connect(self.evt_update_progress)
            # self.worker.update_error_sim_finished.connect(self.evt_update_error_finished)
            self.progressBar.setVisible(True)

            #break

    def runErrorSimulator(self):
        self.inputDNAPath = self.file_path_lineEdit.text()
        while True:
            if not os.path.isfile(self.inputDNAPath):
                print('The chosen input file doesn\'t exist')
                self.msg_box_with_error("The input file you chosen doesn't exist")
                self.file_path_lineEdit.clear()
                break
            else:
                if self.user_defined_copied_checkBox.isChecked():
                    sent_distance_info = self.dist_info
                    if sent_distance_info['type'] == 'continuous' and sent_distance_info['min'] >= sent_distance_info[
                        'max']:
                        self.msg_box_with_error('Min should be < Max in user defined amount of copies')
                        break
                    if sent_distance_info['value'] == '':
                        self.msg_box_with_error('Please enter a value in user defined amount on copies')
                        break
                else:
                    sent_distance_info = None
                self.evyat_path, self.shuffled_path = self.evyat_files_path()
                self.worker = SimulateErrorsWorker(self.general_errors, self.per_base_errors, self.inputDNAPath
                                                   , self.default_radioButton.isChecked(), sent_distance_info, self.evyat_path, self.shuffled_path)
                self.worker.start()
                self.label_progress.setText('Injecting errors, please wait!')
                self.worker.finished.connect(self.evt_worker_finished)
                self.worker.update_progress.connect(self.evt_update_progress)
                # self.worker.update_error_sim_finished.connect(self.evt_update_error_finished)
                self.progressBar.setVisible(True)

                break

    def evt_worker_finished(self):
        self.label_progress.setText('We are done :)')
        self.progressBar.setValue(100)
        self.progressBar.setVisible(False)

    def evt_cluster_worker_finished(self):
        self.evt_worker_finished()

        if platform.system() == "Linux":
            results_f = '/home_nfs/sgamer8/DNAindex' + str(
                self.clustering_index) + '/cluster_output/' + self.chosen_technology + '/' + str(
                self.clustering_index) + '_final_results.txt'
        elif platform.system() == "Windows":
            results_f = 'cluster_output/' + self.chosen_technology + '/' + str(self.clustering_index) + '_final_results.txt'

        text = open(results_f).read()
        self.clustering_results_textEdit.setText(text)


    def evt_update_progress(self, val):
        self.progressBar.setValue(val + 1)

    # def evt_update_error_finished(self, val):
    #     if val == 'error_sim_finished':
    #         self.label_progress.setText('Running reconstruction, please wait!')
    #         self.progressBar.setValue(0)

    def msg_box_with_error(self, error_msg):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(error_msg)
        msg.setWindowTitle("Error!")
        msg.setStandardButtons(QMessageBox.Ok)
        retval = msg.exec_()

    def parse_hist_results(self, working_dir):
        num_clusters = 0
        start_copying = 0
        # x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # y = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
        x = []
        y = []

        source = open(working_dir + '/output.txt', 'r')
        f = open(working_dir + '/histogram.txt', 'w', newline='\n')

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

        source = open(working_dir + '/histogram.txt', 'r')
        for line in source:
            line_list = re.split(r'\t+|\s+', line)  # seperates a string to a list by a delimeter of spaces or tabs
            index = int(line_list[0].strip())
            value = int(line_list[1].strip())
            x.insert(index, index)
            y.insert(index, (value / num_clusters) * 100)
            # y[index] = (value / num_clusters) * 100
        source.close()
        return x, y, num_clusters

    def show_hist_graph_result(self, working_dir):
        import numpy as np
        import matplotlib.pyplot as plt

        # if os.path.isfile('output/histogram.png'):
        #     os.remove('output/histogram.png')

        # self.histogram_img.setText('')

        x, y, num_clusters = self.parse_hist_results(working_dir)

        plt.xticks(x)

        plt.figure()  # creates a new plot, so it doesn't plot a couple in one figure
        plt.scatter(x, y, color='r', zorder=2)
        plt.plot(x, y, color='b', zorder=1)

        plt.title("Reconstruction Errors Histogram")
        plt.xlabel("Number of edit errors")
        plt.ylabel("Fraction of reads")

        plt.savefig(working_dir + '/histogram.png')
        # plt.show()

        # qim = ImageQt('output/histogram.png').copy()
        # pix = QtGui.QPixmap.fromImage(qim)
        # self.histogram_img.setPixmap(pix)
        # self.histogram_img.adjustSize()

        pixmap = QPixmap(working_dir + '/histogram.png')
        self.histogram_img.setPixmap(pixmap)

        # Optional, resize window to image size
        # self.resize(pixmap.width(), pixmap.height())

    def dataReady(self):
        x = str(self.process.readAll(), 'utf-8')
        res = re.split('\r\n', x.strip())
        for i in res:
            if i.isnumeric():
                self.progressBar.setValue(int(i))

    def reconstruction_finished(self):
        if platform.system() == "Linux":
            working_dir = '/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology
        elif platform.system() == "Windows":
            working_dir = 'files/' + self.chosen_technology
        self.label_progress.setText('We are done :)')
        self.progressBar.setVisible(False)
        text = open(working_dir + '/output.txt').read()
        self.reconstruction_output_textEdit.setText(text)
        self.show_hist_graph_result(working_dir)

    def call_reconstruction_alg(self, alg_file_name, working_dir):
        process_path = None
        if platform.system() == "Linux" or platform == "linux2":
            # linux
            process_path = '../reconstruction_algs/' + alg_file_name + "Linux"
        elif platform.system() == "Darwin":
            # OS X
            process_path = '../reconstruction_algs/' + alg_file_name + "Mac"
        elif platform.system() == "Windows":
            # Windows...
            process_path = 'reconstruction_algs/' + alg_file_name + ".exe"
        self.progressBar.setVisible(True)
        self.label_progress.setText('Running reconstruction, please wait!')
        self.process = QProcess(self)
        self.process.setWorkingDirectory(working_dir)
        self.process.start(process_path)
        self.process.readyRead.connect(self.dataReady)
        self.process.finished.connect(self.reconstruction_finished)

    def run_reconstruction_algo(self):
        if platform.system() == "Linux":
            evyat_path = '/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology + '/' + 'evyat.txt'
            working_dir = '/home_nfs/sgamer8/DNAindex' + str(self.clustering_index) + '/files/' + self.chosen_technology
        elif platform.system() == "Windows":
            evyat_path = 'files/' + self.chosen_technology + '/' + 'evyat.txt'
            working_dir = 'files/' + self.chosen_technology

        if not os.path.isfile(evyat_path):
            self.msg_box_with_error('Please run clustering first, or provide the evyat.txt input file')
            self.label_progress.setText('')
            return

        self.label_progress.setText('Running reconstruction, please wait!')
        if self.reconstruction_algo == 'Hybrid Reconstruction Algorithm':
            self.call_reconstruction_alg('Hybrid', working_dir)
        elif self.reconstruction_algo == 'Divider BMA Reconstruction Algorithm':
            self.call_reconstruction_alg('DivBMA', working_dir)
        elif self.reconstruction_algo == 'BMA Look Ahead Reconstruction Algorithm':
            self.call_reconstruction_alg('BMALookahead', working_dir)
        elif self.reconstruction_algo == 'Iterative Reconstruction Algorithm':
            self.call_reconstruction_alg('Iterative', working_dir)
        elif self.reconstruction_algo == 'Mock Reconstruction Algorithm':
            self.call_reconstruction_alg('mockReconstruction', working_dir)
        else:
            self.msg_box_with_error('Please choose a reconstruction algorithm')
            self.label_progress.setText('')
            return

        if not os.path.isfile('output/mock.txt'):
            self.msg_box_with_error('Reconstruction doesn\'t have an output. Try running it again')
            return


class SimulateErrorsWorker(QThread):
    update_progress = pyqtSignal(int)

    # update_error_sim_finished = pyqtSignal(str)

    def __init__(self, general_errors, per_base_errors, input_dna_path, stutter_chosen, dist_info, evyat_path, shuffled_path):
        super(SimulateErrorsWorker, self).__init__()
        self.general_errors = general_errors
        self.per_base_errors = per_base_errors
        self.inputDNAPath = input_dna_path
        self.stutter_chosen = stutter_chosen
        self.dist_info = dist_info
        self.evyat_path = evyat_path
        self.shuffled_path = shuffled_path

    def report_func(self, total_lines, curr_line):
        percent = int(curr_line * 100 // total_lines)
        self.update_progress.emit(percent)

    def run(self):
        # parse value to list if type is vector:
        if self.dist_info is not None and self.dist_info['type'] == 'vector':
            self.dist_info['value'] = ast.literal_eval(self.dist_info['value'])
        error_sim = Simulator(self.general_errors, self.per_base_errors, self.inputDNAPath
                              , self.stutter_chosen, self.dist_info)
        error_sim.simulate_errors(self.report_func, self.evyat_path, self.shuffled_path)
        # self.update_error_sim_finished.emit


class ClusteringWorker(QThread):
    update_progress = pyqtSignal(int)

    def __init__(self, cluster_tech, cluster_index):
        super(ClusteringWorker, self).__init__()
        self.cluster_tech = cluster_tech
        self.cluster_index = cluster_index

    def report_func(self, total_lines, curr_line):
        percent = int(curr_line * 100 // total_lines)
        self.update_progress.emit(percent)

    def run(self):
        start_time = time.time()
        clusters = Clustering(self.cluster_index, self.cluster_tech)
        clusters.create_all_keys()
        clusters.fill_dict_from_shuffled()
        clusters.full_edit_distance(self.report_func)
        clusters.create_evyat_dict()
        [num_of_errors, num_of_false_negative] = clusters.compare_evyat_with_clustering()

        if platform.system() == "Linux":
            results_f = '/home_nfs/sgamer8/DNAindex' + str(self.cluster_index) + '/cluster_output/' + self.cluster_tech + '/' + str(
                self.cluster_index) + '_final_results.txt'
        elif platform.system() == "Windows":
            results_f = 'cluster_output/' + self.cluster_tech + '/' + str(self.cluster_index) + '_final_results.txt'

        with open(results_f, 'w') as results_file:
            print("Run time: %s seconds" % round((time.time() - start_time),2), file=results_file)
            print('Number of false positives: ' + num_of_errors, file=results_file)
            print('Number of false negatives: ' + num_of_false_negative, file=results_file)
            print('Number of thrown strands: ' + str(clusters.total_amount_of_thrown), file=results_file)


if __name__ == '__main__':
    # Create the Qt Application
    app = QApplication(sys.argv)
    # Create and show the form
    dnaSimulator = dnaSimulator()
    dnaSimulator.show()
    # Run the main Qt loop
    sys.exit(app.exec_())
