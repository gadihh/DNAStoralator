import sys
import time

ERRORS_FILE_PATH = 'temp/errors_logger.txt'

# Singelton object.
error_logger = None


def get_logger():
    '''

    :return: ErrorLogger
    '''
    global error_logger

    if not error_logger:
        error_logger = ErrorLogger()

    return error_logger


class ErrorLogger():
    def __init__(self):
        # Clear any content in the file.
        open(ERRORS_FILE_PATH, 'w').close()

    def log_error(self, error_string):
        # Write error to file.
        f = open(ERRORS_FILE_PATH, 'a+')
        f.write(error_string)

        # Write error to sys error output.
        print("There is a problem with the design input either supply a valid csv file or an iupac string.",
              file=sys.stderr)

        f.close()
        time.sleep(0.1)
