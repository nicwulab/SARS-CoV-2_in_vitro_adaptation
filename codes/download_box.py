import requests
from pathlib import Path
import logging

from tqdm import tqdm
import defopt
logging.basicConfig(level=logging.INFO)


class Box:
    def __init__(self, token=''):
        self.baseurl = 'https://api.box.com/2.0'
        self.headers = {'Authorization': 'Bearer %s' %token}

    def list_folder(self, folder_id):
        '''
        API reference: https://developer.box.com/reference/get-folders-id-items/
        '''
        folder_id = int(folder_id)
        r = requests.get(self.baseurl + '/folders/%i/items/' %folder_id, 
                headers = self.headers)
        if r.status_code == 200:
            return r.json()['entries']
        else:
            raise Exception('Return code: {}'.format(r.status_code))

    def download_file(self, file_id, out_file):
        '''
        API reference: https://developer.box.com/reference/get-files-id-content/

        download function: https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests
        '''

        file_id = int(file_id)
        url = self.baseurl + '/files/%i/content/' %file_id
        # Streaming, so we can iterate over the response.
        response = requests.get(url, stream=True, headers = self.headers)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        block_size = 1024 #1 Kibibyte
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc=str(out_file))
        mode = 'wb' 
        with open(out_file, mode) as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            raise Exception("ERROR, something went wrong")


def main(token: str, outpath: str, folder_id: int = 129452943094, filename: str='trimmed.bam'):
    '''
    Download bam files from Box, folder structure:
    - folder/sample1/trimmed.bam
    - folder/sample2/trimmed.bam
    - folder/sample3/trimmed.bam

    The downloaded files will be:
    - outpath/sample1/timmed.bam
    - outpath/sample2/timmed.bam
    - outpath/sample3/timmed.bam

    The folder is defined by folder_id

    :param token: developer token from box custom app 
    :param outpath: where to store downloaded files
    :param folder_id:  Box folder id 
    :param filename:  filename from the sub-folder to be downloaded

    '''
    box = Box(token=token)
    outpath = Path(outpath)
    for sample in box.list_folder(folder_id):
        if sample['type'] == 'folder':
            for sample_file in box.list_folder(sample['id']):
                if sample_file['name'] == filename:
                    out_folder = outpath / sample['name']
                    out_folder.mkdir(exist_ok=True)
                    out_file_name = out_folder / filename
                    box.download_file(sample_file['id'], out_file_name)
                    logging.info('Downloaded %s' %out_file_name)

if __name__ == '__main__':
    defopt.run(main)
