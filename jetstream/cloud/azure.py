import os
import json
import subprocess
from datetime import datetime

from .base import CloudStorageSession


class AzureStorageSession(CloudStorageSession):
    def __init__(self, az_sif_path, az_storage_account_name, az_storage_account_key='',
                 create_temp_container=True, keep_temp_container=False):
        super().__init__()
        self.az_sif_path = az_sif_path
        self.keep_temp_container = keep_temp_container
        
        # Log into azure account
        # print('Logging into Azure CLI')
        # self.provider_login(az_username, az_password)
        
        # Get an account key to be used in subsequent transactions
        self.storage_account_name = az_storage_account_name
        # print('Getting Account Key')
        self.storage_account_key = az_storage_account_key
        # print(f'Account key: {self.account_key}')
        
        if create_temp_container:
            print('Creating Temp Container')
            self.create_temp_container()
        
    # def _account_key(self):
    #     cmd = (f"""
    #         singularity exec {self.az_sif_path} az storage account keys list 
    #             -g {self.resource_group} 
    #             -n {self.storage_account_name}
    #     """).split()
    #     return json.loads(subprocess.check_output(cmd).decode())[0]['value']
    
    
    
    
    def get_account_key(self, az_username, az_password, resource_group, storage_account_name):
        login_cmd = (f"""
            singularity exec {self.az_sif_path} az login -u {az_username} -p "{az_password}"
        """).split()
        subprocess.call(login_cmd)
        
        get_key_cmd = (f"""
            singularity exec {self.az_sif_path} az storage account keys list 
                -g {self.resource_group} 
                -n {self.storage_account_name}
        """).split()
        self.account_key = json.loads(subprocess.check_output(get_key_cmd).decode())[0]['value']
        return self.account_key
        
    
    # def provider_login(self, az_username, az_password):
    #     login_cmd = (f"""
    #         singularity exec {self.az_sif_path} az login -u {az_username} -p "{az_password}"
    #     """).split()
    #     subprocess.call(login_cmd)
    
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
                singularity exec {self.az_sif_path} az storage container create 
                    -n {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            self._container_created = response['created']
    
    async def _blob_container_interact(self, interaction, filepath, blobpath=None, container=None):
        container = container or self._temp_container_name
        
        # If blob path on the container isn't given, make it the same as the filepath
        blobpath = blobpath or os.path.basename(filepath)
        
        print(f'{interaction.capitalize()}ing file {filepath}')
        
        # TODO that --bind /mnt:/mnt is a hack for now and needs to be revisited later
        cmd = (f"""
            singularity exec --bind /mnt:/mnt {self.az_sif_path} az storage blob {interaction} 
                --container-name {container} 
                --file {filepath} 
                --name {blobpath} 
                --account-name {self.storage_account_name} 
                --account-key {self.storage_account_key} 
        """).split()
        # print(f'called: {cmd}')
        # await self.subprocess_sh(cmd)
        await self.subprocess_sh(' '.join(cmd))
        # print('finished')
        # subprocess.check_output(cmd)
    
    async def upload_blob(self, filepath, blobpath=None, container=None, force=False):
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        blobpath = blobpath or os.path.basename(filepath)
        if not force:
            cmd = (f"""
                singularity exec {self.az_sif_path} az storage blob exists 
                    --container-name {self._temp_container_name}
                    --name {blobpath} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if response.get('exists'):
                print(f'Blob {blobpath} already exists')
                return
            
        await self._blob_container_interact(
            interaction='upload',
            filepath=filepath,
            blobpath=blobpath,
            container=container
        )
        
    async def download_blob(self, filepath, blobpath=None, container=None):
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        await self._blob_container_interact(
            interaction='download',
            filepath=filepath,
            blobpath=blobpath,
            container=container
        )
    
    def download_blobs_as_bash(self, blobs, output_paths, container=None):
        print('download blobs as bash')
        cmd = ''
        for blob, output_path in zip(blobs, output_paths):
            cmd += (f'az storage blob download --container-name {container} --file {output_path} --name {blob}'
                     '--account-name {self.storage_account_name} --account-key {self.storage_account_key};')
        return cmd
    
    def upload_blobs_as_bash(self):
        pass
            
        
    def close(self):
        if not self.storage_account_key:
            self._no_storage_account_key()
        
        if not self.keep_temp_container and self._container_created:
            cmd = (f"""
                singularity exec {self.az_sif_path} az storage container delete 
                    --name {self._temp_container_name} 
                    --account-name {self.storage_account_name} 
                    --account-key {self.storage_account_key}
            """).split()
            response = json.loads(subprocess.check_output(cmd).decode())
            if not response.get('deleted'):
                print('The container could NOT be deleted, do it manually to avoid useage charges')