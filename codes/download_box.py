import requests
from tqdm import tqdm
import defopt

TOKEN='koDT1eUsVlgzymO6ZwqMh7ARpx4FfrKN'



class Box:
    def __init__(self, token=''):
        self.TOLEN = token
        self.baseurl = 'https://api.box.com/2.0'
        self.headers = {'Authorization': 'Bearer %s' %TOKEN}

    def list_folder(self, folder_id):
        '''
        API reference: https://developer.box.com/reference/get-folders-id-items/
        '''
        folder_id = int(folder_id)
        r = requests.get(self.baseurl + '/folders/%i/items/' %folder_id, 
                headers = self.headers)
        return r.json()['entries']

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
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc=out_file)
        with open(out_file, 'wb') as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            print("ERROR, something went wrong")

def main(token: str, outpath: str, folder_id: int = 129452943094):
    '''
    Download bam files from Box, folder structure:
    - folder/sample1/trimmed.bam
    - folder/sample2/trimmed.bam
    - folder/sample3/trimmed.bam

    The downloaded files will be:
    - outpath/sample1.timmed.bam
    - outpath/sample2.timmed.bam
    - outpath/sample3.timmed.bam

    The folder is defined by folder_id

    :param token: developer token from box custom app 
    :param outpath: where to store downloaded files
    :param folder id:  Box folder id 

    '''
    box = Box(token=token)
    for sample in box.list_folder():
        if sample['type'] == 'folder':
            for sample_file in box.list_folder(sample['id']):
                if sample_file['name'] == 'trimmed.bam':
                    bam_name = outpath + '/' + sample['name'] + '.trimmed.bam'
                    box.download_file(sample_file['id'], bam_name)
                


if __name__ == '__main__':
    defopt.run(main)
