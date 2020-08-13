import os
import json
import subprocess
from datetime import datetime
import urllib

from .base import CloudStorageSession, is_remote_uri, path_conversion


class AzureStorageSession(CloudStorageSession):
    """
    Provider for Azure blob storage.
    """
    config_key = 'azure_params'
    
    def __init__(self, az_storage_account_name, az_storage_account_key='', az_sif_path=None, 
                 create_temp_container=True, keep_temp_container=False, blob_container=None):
        """
        The primary method for accessing the ``az`` binary is through a singularity container, but if no path to a 
        container is given, it is assumed that ``az`` is in PATH. The user must also provide a storage account name and 
        the associated account key.
        
        :param az_storage_account_name: str A valid Azure storage account name
        :param az_storage_account_key: str The corresponding Azure storage account key
        :param az_sif_path: str If not provided, will download the official latest container from Microsoft
        :param create_temp_container: bool Whether to create a new container for this run
        :param keep_temp_container: bool Whether to keep the temporary container after the run completes
        # TODO The above parameter currenly has no effect, because I don't know where in the jetstream lifecycle 
          I can call the AzureStorageSession.close() function, which would clean up the temp container
        :param blob_container: str If given, this Azure storage container will be used instead of creating a new 
        temp container
        """
        super().__init__(container=blob_container)
        self.az_sif_path = az_sif_path or 'docker://mcr.microsoft.com/azure-cli'
        self.keep_temp_container = keep_temp_container
        singularity_available = subprocess.call(['which', 'singularity']) == 0
        self.az_execution_call = f'singularity exec --bind /mnt:/mnt {self.az_sif_path} az' if singularity_available else 'az'
        
        # Get an account key to be used in subsequent transactions
        self.storage_account_name = az_storage_account_name
        self.storage_account_key = az_storage_account_key
        
        if create_temp_container:
            print('Creating Temp Container')
            self.create_temp_container()
    
    def remote_download_cmd(self, remote_downloads_blobs=None):
        remote_download_cmd_template = (
            'if [[ ! -f "{blobpath}" ]];then mkdir -p {blobpath_dirname}; '
            'singularity exec {az_sif_path} az storage blob download --name {blobpath} --file {blobpath} '
            '--container-name {container_name} --account-name {account_name} --account-key {account_key};fi\n'
        )
        remote_download_cmd = ''
        for remote_download_blobpath in remote_downloads_blobs or list():
            if is_remote_uri(remote_download_blobpath):
                # TODO There is currently no mechanism to tell wget where to put the file
                url_filepath = urllib.parse.urlparse(remote_download_blobpath).path
                remote_download_cmd += 'if [[ ! -f "{}" ]];then wget {};fi\n'.format(
                    os.path.basename(url_filepath),
                    remote_download_blobpath
                )
            else:
                remote_download_blobpath = path_conversion(remote_download_blobpath)
                remote_download_cmd += remote_download_cmd_template.format(
                    blobpath=remote_download_blobpath,
                    blobpath_dirname=os.path.dirname(remote_download_blobpath),
                    az_sif_path=self.az_sif_path,
                    container_name=self.container,
                    account_name=self.storage_account_name,
                    account_key=self.storage_account_key
                )
        
        return remote_download_cmd
    
    def remote_upload_cmd(self, remote_uploads_filepaths=None):
        remote_upload_cmd_template = (
            'singularity exec {az_sif_path} az storage blob upload --name {blobpath} --file {blobpath} '
            '--container-name {container_name} --account-name {account_name} --account-key {account_key};\n'
        )
        remote_upload_cmd = ''
        for remote_upload_filepath in remote_uploads_filepaths or list():
            remote_upload_filepath = path_conversion(remote_upload_filepath)
            remote_upload_cmd += remote_upload_cmd_template.format(
                blobpath=remote_upload_filepath,
                blobpath_dirname=os.path.dirname(remote_upload_filepath),
                az_sif_path=self.az_sif_path,
                container_name=self.container,
                account_name=self.storage_account_name,
                account_key=self.storage_account_key
            )
        return remote_upload_cmd
    
    def get_account_key(self, az_username, az_password, resource_group, storage_account_name):
        login_cmd = (f"""
            {self.az_execution_call} login -u {az_username} -p "{az_password}"
        """).split()
        subprocess.call(login_cmd)
        
        get_key_cmd = (f"""
            {self.az_execution_call} storage account keys list 
                -g {self.resource_group} 
                -n {self.storage_account_name}
        """).split()
        self.account_key = json.loads(subprocess.check_output(get_key_cmd).decode())[0]['value']
        return self.account_key
    
    def _no_storage_account_key(self):
        raise ValueError(
            'Azure storage account key is not set, either provide it during instantiation or call '
            'AzureStorageSession.get_account_key()'
        )
    
    def create_temp_container(self):
        """
        self._temp_container_name comes from the superclass.
        """
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        if not self._container_created:
            # Container needs to be created
            cmd = (f"""
                {self.az_execution_call} storage container create 
                    -n {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            self._container_created = response['created']
    
    def upload_blob(self, filepath, blobpath=None, container=None, force=False):
        if not self.storage_account_key:
            self._no_storage_account_key()
        

        if is_remote_uri(filepath):
            print(f'Blob {blobpath} is a remote URI, nothing to upload')
            return
        
        blobpath = blobpath or os.path.basename(filepath)
        if not force:
            cmd = (f"""
                {self.az_execution_call} storage blob exists 
                    --container-name {self._temp_container_name} 
                    --name {blobpath} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if response.get('exists'):
                print(f'Blob {blobpath} already exists')
                return
        
        # Do the upload stuffs
        container = container or self._temp_container_name
        
        # If blob path on the container isn't given, make it the same as the filepath
        blobpath = blobpath or os.path.basename(filepath)
        
        print(f'uploading file {filepath}')
        
        # TODO that --bind /mnt:/mnt is a hack for now and needs to be revisited later
        cmd = (f"""
            {self.az_execution_call} storage blob upload
                --container-name {container} 
                --file {filepath} 
                --name {blobpath} 
                --account-name {self.storage_account_name} 
                --account-key {self.storage_account_key} 
        """).split()
        subprocess.check_output(cmd)
    
        
    def download_blob(self, filepath, blobpath=None, container=None):
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        # Do the download things
        container = container or self._temp_container_name
        
        # If blob path on the container isn't given, make it the same as the filepath
        blobpath = blobpath or os.path.basename(filepath)
        
        print(f'downloading file {filepath}')
        
        # TODO that --bind /mnt:/mnt is a hack for now and needs to be revisited later
        cmd = (f"""
            {self.az_execution_call} storage blob download
                --container-name {container} 
                --file {filepath} 
                --name {blobpath} 
                --account-name {self.storage_account_name} 
                --account-key {self.storage_account_key} 
        """).split()
        subprocess.check_output(cmd)
        
    def close(self):
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        if not self.keep_temp_container and self._container_created:
            cmd = (f"""
                {self.az_execution_call} storage container delete 
                    --name {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if not response.get('deleted'):
                print('The container could NOT be deleted, do it manually to avoid useage charges')
