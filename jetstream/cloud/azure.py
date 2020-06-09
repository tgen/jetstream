import os
import json
import subprocess
from datetime import datetime

# from .base import CloudStorageSession


class AzureStorageSession(CloudStorageSession):
    def __init__(self, az_sif_path, az_username, az_password, resource_group, storage_account_name,
                 create_temp_container=True, keep_temp_container=False):
        super().__init__()
        self.az_sif_path = az_sif_path
        self.keep_temp_container = keep_temp_container
        
        # Log into azure account
        print('Logging into Azure CLI')
        self.provider_login(az_username, az_password)
        
        # Get an account key to be used in subsequent transactions
        self.resource_group = resource_group
        self.storage_account_name = storage_account_name
        print('Getting Account Key')
        self.account_key = self._account_key()
        print(f'Account key: {self.account_key}')
        
        if create_temp_container:
            print('Creating Temp Container')
            self.create_temp_container()
        
    def _account_key(self):
        cmd = (f"""
            singularity exec {self.az_sif_path} az storage account keys list 
                -g {self.resource_group} 
                -n {self.storage_account_name}
        """).split()
        return json.loads(subprocess.check_output(cmd).decode())[0]['value']
    
    def provider_login(self, az_username, az_password):
        login_cmd = (f"""
            singularity exec {self.az_sif_path} az login -u {az_username} -p "{az_password}"
        """).split()
        subprocess.call(login_cmd)
    
    def create_temp_container(self):
        """
        self._temp_container_name comes from the superclass.
        """
        if not self._container_created:
            # Container needs to be created
            cmd = (f"""
                singularity exec {self.az_sif_path} az storage container create 
                    -n {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            self._container_created = response['created']
    
    def _blob_container_interact(self, interaction, filepath, blobpath=None, container=None):
        container = container or self._temp_container_name
        
        # If blob path on the container isn't given, make it the same as the filepath
        blobpath = blobpath or filepath
        
        cmd = (f"""
            singularity exec {self.az_sif_path} az storage blob {interaction} 
                --container-name {container} 
                --file {filepath} 
                --name {blobpath} 
                --account-name {self.storage_account_name} 
                --account-key {self.account_key}
        """).split()
        subprocess.check_output(cmd)
    
    def upload_blob(self, filepath, blobpath=None, container=None, force=False):
        blobpath = blobpath or filepath
        if not force:
            cmd = (f"""
                singularity exec {self.az_sif_path} az storage blob exists 
                    --container-name {self._temp_container_name}
                    --name {blobpath} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if response.get('exists'):
                print(f'Blob {blobpath} already exists')
                return
            
        self._blob_container_interact(
            interaction='upload',
            filepath=filepath,
            blobpath=blobpath,
            container=container
        )
        
    def download_blob(self, filepath, blobpath=None, container=None):
        self._blob_container_interact(
            interaction='download',
            filepath=filepath,
            blobpath=blobpath,
            container=container
        )
        
    def close(self):
        if not self.keep_temp_container and self._container_created:
            cmd = (f"""
                singularity exec {self.az_sif_path} az storage container delete 
                    --name {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if not response.get('deleted'):
                print('The container could NOT be deleted, do it manually to avoid useage charges')