from datetime import datetime
from abc import ABC, abstractmethod

class CloudStorageSession(ABC):
    def __init__(self):
        self._temp_container_name = 'jetstream-temp-{unique_datestamp}'.format(
            unique_datestamp=datetime.now().strftime('%Y%m%d%H%M%S%f')
        )
        self._container_created = False
    
    # @abstractmethod
    # def provider_login(self, *args):
    #     raise NotImplementedError
    
    @abstractmethod
    def upload_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError
    
    @abstractmethod
    def download_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError