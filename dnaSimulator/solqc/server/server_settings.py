'''

this file contains all configurations and settings for server
'''
import os


class Config:

    def __init__(self):
        self.NONE = "none"
        self.ERROR_CONSTANT = 'error'
        self.DOWNLOAD_ERROR = 'download error'
        self.DELETE_FILES = True
        self.SERVER_FOLDER = 'server/'
        # not implemented yet but maybe is necessary
        self.DELETE_REPORTS = True
        self.cwd = os.getcwd()
        self.starting_path = os.getcwd()
        self.DEFAULT_DATA_PATH = os.path.join(self.cwd,  "data")
        self.DELIVERABLE_PATH = 'deliverable/'
        self.TEMP_PATH = os.path.join(self.cwd,  "temp")
        self.prefix = "input_data"
        self.UPLOAD_FOLDER = self.SERVER_FOLDER + "uploaded_files/"
        self.READS_FROM_ZIP =self.UPLOAD_FOLDER + "temp_reads_directory/"
        self.MAX_FILE_SIZE = 20*1024*1024*1024  # 30 GB

    def create_directories(self):
        # create upload file folder if it doesnt exists
        if not os.path.isdir(self.UPLOAD_FOLDER):
            os.makedirs(self.UPLOAD_FOLDER)

        # create server_data folder if it doesnt exists
        if not os.path.isdir(self.DEFAULT_DATA_PATH):
            os.makedirs(self.DEFAULT_DATA_PATH)

        # create temp folder to unzip files to
        if not os.path.isdir(self.READS_FROM_ZIP):
            os.makedirs(self.READS_FROM_ZIP)

    def get_data_path(self):
        return self.DEFAULT_DATA_PATH

    def get_new_data_path(self):
        return self.DEFAULT_DATA_PATH + "/random"

    def set_data_path(self,data_path):
        self.DEFAULT_DATA_PATH = data_path


config = Config()


