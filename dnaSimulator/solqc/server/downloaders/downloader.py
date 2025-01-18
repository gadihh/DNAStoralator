import wget
import os
from enum import Enum
from server.downloaders import google_downloader


class Provider(Enum):
    DEFAULT = 1
    GOOGLE_DRIVE = 2

GOOGLE = 'Google'

#######################################################################################################################
#                                                           API
#######################################################################################################################
def download_by_link(link, destination):
    """
    Downloads a file by the given link and stores it in the given destination.
    :param link: The link to downloadble file.
    :param destination: A string of the like - "/path/filename.fastq"
    :return: True if all went well, False otherwise.
    """
    provider = resolve_provider(link)
    if provider == Provider.GOOGLE_DRIVE:
        print("Using Google drive to download")

        success = google_downloader.download(link, destination)

        return success

    elif provider == Provider.DEFAULT:
        print("Using wget to download")
        # ret = os.system("wget  https://www.dropbox.com/s/xfbl0lccd614m8a/Archive.zip?dl=0")
        ret = os.system(f'wget -O {destination} {link}')
        if ret == 0:
            return True

        return False

    # shoud not get to here
    return False
#######################################################################################################################
#                                                        Private
#######################################################################################################################

def resolve_provider(link):

    # TODO change to a better regex
    if("google" in link):
        return Provider.GOOGLE_DRIVE

    return Provider.DEFAULT

#######################################################################################################################
#                                                        Test
#######################################################################################################################
if __name__ == '__main__':
    print(download_by_link("https://www.dropbox.com/s/xfbl0lccd614m8a/Archive2.zip", "test.zip"))