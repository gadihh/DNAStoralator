import requests
import re

#######################################################################################################################
#                                                           API
#######################################################################################################################

def download(link, destination):
    """
    Download a file from google drive using the sharable link.
    The file will be downloaded and saved in the destination.
    Destanation should include the name of the file as well as the path.
    :param link: Google drive generated sharable link
    :param destination: A string of the like - "/path/filename.fastq"
    :return: True if all went well, False otherwise.
    """
    file_id = extract_id_from_link(link)

    if not file_id:
        exit_with_error()
        return False
    else :
        download_file_from_google_drive(file_id, destination)
        return True

#######################################################################################################################
#                                                         Private
#######################################################################################################################
def exit_with_error():
    print("Failed to download file from google drive")


def extract_id_from_link(sharable_link):
    id_pattern = r'(/[-\w]{25,}/)'

    result = re.search(id_pattern, sharable_link)
    if  result: # first possible sharable link
    # Remove the starting and trailing dashes(/)
        return result.group()[1:-1]
        # https://drive.google.com/open?id=12s3XC8qKDyQPN0josFp1QeJLkerYjP9Q

    id_pattern = r'\?id=(.*)'
    result = re.search(id_pattern, sharable_link)
    if result:# second possibility for sharable link
        return result.group(1)

    # in case of no matches
    print("Invalid sharable link to google drive")
    return False


def download_file_from_google_drive(id, destination):
    def get_confirm_token(response):
        for key, value in response.cookies.items():
            if key.startswith('download_warning'):
                return value

        return None

    def save_response_content(response, destination):
        CHUNK_SIZE = 32768

        with open(destination, "wb") as f:
            for chunk in response.iter_content(CHUNK_SIZE):
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)

    URL = "https://docs.google.com/uc?export=download"

    session = requests.Session()

    response = session.get(URL, params = { 'id' : id }, stream = True)
    token = get_confirm_token(response)

    if token:
        params = { 'id' : id, 'confirm' : token }
        response = session.get(URL, params = params, stream = True)

    save_response_content(response, destination)


#######################################################################################################################
#                                                         Test
#######################################################################################################################
if __name__ == "__main__":
    # Extract id from link
    link = "https://drive.google.com/file/d/1isl9DrBVmex06MwWSvHhHTCrrmP0wcXq/view?usp=sharing"
    link2 = "https://drive.google.com/open?id=1oDEinguMQrZdW2oG-IzZQVc1gGvnGe3P"
    link3 = "https://drive.google.com/open?id=12s3XC8qKDyQPN0josFp1QeJLkerYjP9Q"
    file_id = extract_id_from_link(link3)

    # Destination for downloaded file
    destination = "photo.jpg"

    print(download_file_from_google_drive(file_id, destination))
