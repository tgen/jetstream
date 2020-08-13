from datetime import datetime
from abc import ABC, abstractmethod
import asyncio
from asyncio import create_subprocess_shell
import os
import urllib
import time

CLOUD_STORAGE_PROVIDERS = {
    'azure': 'jetstream.cloud.azure.AzureStorageSession'
}


def is_remote_uri(uri):
    return urllib.parse.urlparse(uri).scheme != ''


def dummy_dot_in_path(path):
    return '.' in path.split(os.path.sep)
        

def path_conversion(path):
    if is_remote_uri(path):
        return path
    
    if dummy_dot_in_path(path):
        # local: /absolute/path/./to/file ----> remote: current/dir/to/file
        # local: relative/path/./to/file ----> remote: current/dir/to/file
        # local: current/dir/to/file <---- remote: /absolute/path/./to/file
        # local: current/dir/to/file <---- remote: relative/path/./to/file
        return path.rsplit(f'.{os.path.sep}')[-1]
    
    if os.path.isabs(path):
        # local: /absolute/path/to/file ----> remote: /absolute/path/to/file
        # local: /absolute/path/to/file <---- remote: /absolute/path/to/file
        return path
    
    # local: relative/path/to/file ----> remote: current/dir/file
    # local: current/dir/file <---- remote: relative/path/to/file
    return path.split(os.path.sep)[-1]


def blob_inputs_to_remote(inputs, cloud_storage, blob_basename=False):
    elapsed_times = list()
    for input_file in inputs:
        try:
            start = time.time()
            cloud_storage.upload_blob(
                input_file,
                blobpath=os.path.basename(input_file) if blob_basename else path_conversion(input_file)
            )
                
            elapsed_times.append({
                'name': input_file,
                'size': os.stat(input_file).st_size,
                'time': time.time() - start
            })
        except:
            elapsed_times.append({'name': input_file, 'size': -1, 'time': -1})
    return elapsed_times
    

def blob_outputs_to_local(outputs, cloud_storage, blob_basename=False):
    elapsed_times = list()
    for output_file in outputs:
        try:
            start = time.time()
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            cloud_storage.download_blob(
                output_file,
                blobpath=os.path.basename(output_file) if blob_basename else path_conversion(output_file)
            )
            elapsed_times.append({
                'name': output_file,
                'size': os.stat(output_file).st_size,
                'time': time.time() - start
            })
        except Exception as e:
            elapsed_times.append({'name': output_file, 'size': -1, 'time': -1})
    return elapsed_times


class CloudStorageSession(ABC):
    """
    Base for all cloud storage provider specifications.
    """
    def __init__(self, container=None):
        self.container = container or 'jetstream-temp-{unique_datestamp}'.format(
            unique_datestamp=datetime.now().strftime('%Y%m%d%H%M%S%f')
        )
        self._temp_container_name = self.container  # Backward compatibility
        self._container_created = False
    
    @abstractmethod
    def upload_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError
    
    @abstractmethod
    def download_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError
