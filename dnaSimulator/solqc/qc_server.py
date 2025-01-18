from flask import Flask, render_template,session, request,jsonify, redirect, url_for,send_from_directory, send_file, make_response
import threading
import webbrowser
from server.files_manager import *
from server.server_settings import *
from main import start_analyzer
import os, sys

app = Flask(__name__, root_path=config.SERVER_FOLDER, static_folder=os.path.join(os.getcwd(), config.SERVER_FOLDER,'static'))
app.secret_key = os.urandom(30)

def run_analyzer(app,params):
    start_analyzer(params)
    return "done"



@app.route('/')
def home_page():
    if  'current_data_path' not in session:
        session['current_data_path'] = config.DEFAULT_DATA_PATH

    design_files = get_files_with_extension(".csv", session['current_data_path'])
    fastq_files = get_files_with_extension(".fastq",session['current_data_path'])
    config_files = get_files_with_extension(".json",session['current_data_path'])

    # load configuration file from json
    with open(config.SERVER_FOLDER +"server_config_json.json", "r") as read_file:
        server_config = json.load(read_file)

    return render_template('upload_home_page.html', str=str, fastq_files=fastq_files, design_files=design_files,
                           config_files=config_files, server_config=server_config, data_path = session['current_data_path'])



@app.route('/processing/<time_stamp>')
def analyzing_completed(time_stamp):
    # go through all threads and see if thread with same time_stamp exists
    # i.e check if we are still working on generating pdf file
    for mthread in threading.enumerate():
        if mthread.getName() == time_stamp:
            return "still working"
    #     else

    if config.DELETE_FILES:
        delete_all_files_in_directory(config.UPLOAD_FOLDER)
        delete_all_files_in_directory(config.READS_FROM_ZIP)

    # if the thread terminated check if .pdf file is generated
    date = datetime.utcfromtimestamp(int(time_stamp)/1000).strftime('%Y-%m-%d_%H-%M-%S')
    reports_path = os.path.join(config.DELIVERABLE_PATH, date)
    reports = []
    for file in os.listdir(reports_path):
        reports.append(file)

    error_path = os.path.join(config.TEMP_PATH, "errors_logger.txt")
    if os.path.isfile(error_path):
        with open(error_path, 'r') as error_file:
            error_log = error_file.read()

        delete_all_files_in_directory(config.TEMP_PATH)

        return error_log

    if reports:
        return render_template('results.html',date=date,reports=reports)

    return config.ERROR_CONSTANT


def report_page(reports):
    render_template('results.html', reports= reports)


@app.route('/deliverable/<path:filename>',  methods=['GET', 'POST'])
def deliverable_static(filename):
    path = os.path.join(config.cwd, config.DELIVERABLE_PATH)
    return send_from_directory(path, filename)


@app.route('/processing')
def processing_file():
    return render_template('processing_file.html')


@app.route('/uploading', methods=['GET', 'POST'])
def upload_file():

    if request.method == 'POST':
        try:
            time_stamp = request.values.get("time_stamp")
            params = process_data_from_post(request.form, request.files, time_stamp, session['current_data_path'])
            if params == config.DOWNLOAD_ERROR:
                return " we had a problem downloding from the link you supplied, please check the link or upload a file in a different way"
            else:
                thread = threading.Thread(target=run_analyzer, name=time_stamp, args=(app,params))
                thread.start()
                return render_template("processing_file.html", time_stamp=time_stamp)
        except Exception as inst:
            print(inst)
            print('baaad')
            return "we had an error"


if __name__ == '__main__':

    if len(sys.argv) == 2:
        config.set_data_path(sys.argv[1])

    config.create_directories()

    app.config['MAX_CONTENT_LENGTH'] = config.MAX_FILE_SIZE
    app.config['UPLOAD_FOLDER'] = config.DELIVERABLE_PATH
    webbrowser.open('http://localhost:5000/', new=2)
    app.run(host='0.0.0.0', port=5000)
