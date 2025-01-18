import os, sys
import json
from server.downloaders import downloader, google_downloader
from werkzeug.utils import secure_filename
from server.server_settings import *
from datetime import datetime
from pathlib import Path
import zipfile


def unzip_files(path_to_file, destination):
    with zipfile.ZipFile(path_to_file, "r") as zip_ref:
        zip_ref.extractall(destination)

def delete_all_files_in_directory(path):
    for f in os.listdir(path):
        file_path = os.path.join(path, f)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
        except Exception as e:
            print(e)


def get_files_with_extension(extension, dir_path):
    all_files = os.listdir(dir_path)
    return list(filter(lambda f: f.endswith(extension), all_files))


def create_reads_file(list_of_fastq,location, time_stamp):
    read_file_name = str(time_stamp) + 'read_file' + '.txt'
    read_file_path = os.path.join(config.UPLOAD_FOLDER, read_file_name)

    with open(read_file_path, "w") as read_file:
        for f in list_of_fastq:
            if f.startswith(".") or not f.endswith(".fastq"):
                continue
            f_path = os.path.join(location, f)
            read_file.write(f_path + '\n')
    return read_file_path
def save_file(from_server_data, user_upload_data, file_name, time_stamp, current_data_path):
    """

    this functions return path name and saves file if user uploaded

    :param from_server_data:
    :param user_upload_data:
    :param file_name:
    :param time_stamp:
    :return:
    """
    from_server_name = file_name + "_server"
    user_upload_name = file_name + "_upload"
    url_name = file_name + "_url"

    if from_server_data.get(from_server_name):  # if user chose file from server
        file_path = os.path.join(current_data_path, from_server_data[from_server_name])

    if user_upload_data.get(user_upload_name):  # if the user chose to upload file

        file = user_upload_data[user_upload_name]
        file_name = secure_filename(time_stamp + file.filename)
        file_path = os.path.join(config.UPLOAD_FOLDER, file_name)
        file.save(file_path)
    if from_server_data.get(url_name):
        file_path = os.path.join(config.UPLOAD_FOLDER, file_name)
        if not downloader.download_by_link(from_server_data[url_name], file_path):
            return config.DOWNLOAD_ERROR


    return file_path

def process_flags_for_anlyzer(form_data):
    with open(config.SERVER_FOLDER +"server_config_json.json", "r") as read_file:
        server_config = json.load(read_file)
    flags = []
    checked_flags = []
    for inst in server_config:
        flag = server_config[inst]['params']
        options = form_data.getlist(flag)
        if options and flag not in checked_flags:
            checked_flags.append(flag)
            flags.append(flag)
            flags.extend(options)


    return flags


def process_data_from_post(form_data, user_upload_data, time_stamp, current_data_path):

    """
    This method receives post data from the user and does with it ....

    :param form_data:
    :param user_upload_data:
    :param time_stamp:
    :return: params in format to run star_analyzer
    """

    read_file_path = config.NONE
    zip_file_path = config.NONE

    # if user uploads an IUPAC then we don't need to save a design file.
    if form_data.get("IUPAC"):
        design_file_path = form_data['IUPAC']
    else:
        design_file_path = save_file(form_data, user_upload_data, "design_file", time_stamp, current_data_path)

    config_file_path = save_file(form_data, user_upload_data, "config_file", time_stamp, current_data_path)

    fastq_files_server = form_data.getlist('fastq_file_server')
    csv_file_server = form_data.get("csv_file_server")
    if csv_file_server:
        read_file_path = os.path.join(config.DEFAULT_DATA_PATH,csv_file_server)
    if fastq_files_server:  # if user chose  fastq files from server
        read_file_path = create_reads_file(fastq_files_server, current_data_path, time_stamp)

    #     in case a user uploads a zip\fastq file, it unzips it, and create appropriate reads file.
    if user_upload_data.get("zip_file_upload") :
        zip_file = user_upload_data["zip_file_upload"]
        zip_file_name = secure_filename(time_stamp + zip_file.filename)
        zip_file_path = os.path.join(config.UPLOAD_FOLDER, zip_file_name)
        zip_file.save(zip_file_path)
        delete_all_files_in_directory(config.READS_FROM_ZIP)
        # if user uploaded a zip file
        if zip_file_name.endswith(".zip"):
            unzip_files(zip_file_path, config.READS_FROM_ZIP)
            read_file_path = create_reads_file(os.listdir(config.READS_FROM_ZIP), config.READS_FROM_ZIP, time_stamp)

#     if user uploaded a csv file (after matching)
        elif zip_file_name.endswith(".csv"):
            read_file_path = zip_file_path

# if user uploaded a single fastq file
        else:
            read_file_path = create_reads_file([zip_file_name], config.UPLOAD_FOLDER, time_stamp)
    # in the case gives a url to download from, downloads the file, uzip it and creates a reads file.
    if form_data.get("zip_file_url"):
        # if reads is after matching we expect a csv file
        if(form_data.get('after_matching')):
            read_file_path = os.path.join(config.UPLOAD_FOLDER,time_stamp + "reads.csv")
            if not downloader.download_by_link(form_data["zip_file_url"], read_file_path):
                read_file_path = config.DOWNLOAD_ERROR
        #in case reads are before allignment
        else:
            zip_file_path = os.path.join(config.UPLOAD_FOLDER, "fastq_files.zip")
            if downloader.download_by_link(form_data["zip_file_url"], zip_file_path):
                delete_all_files_in_directory(config.READS_FROM_ZIP)
                unzip_files(zip_file_path, config.READS_FROM_ZIP)
                read_file_path = create_reads_file(os.listdir(config.READS_FROM_ZIP), config.READS_FROM_ZIP, time_stamp)
            else:
                read_file_path = config.DOWNLOAD_ERROR


    # change date format
    date = datetime.utcfromtimestamp(int(time_stamp)/1000).strftime('%Y-%m-%d_%H-%M-%S')
    files_args = ['-r', read_file_path, '-d', design_file_path, '-time', date, '-c', config_file_path]
    flag_args = process_flags_for_anlyzer(form_data)
    args = files_args + flag_args

    # check that download worked properly
    if config.DOWNLOAD_ERROR in files_args:
        args = config.DOWNLOAD_ERROR

    return args
