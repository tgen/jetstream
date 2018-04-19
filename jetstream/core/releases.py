import logging
import os
import shutil
import tarfile
import tempfile
from getpass import getpass

import requests
from requests.auth import HTTPBasicAuth

RELEASE_DIR = os.environ.get('JETSTREAM_RELEASES_DIR')
log = logging.getLogger(__name__)


def _basic_auth():
    auth = HTTPBasicAuth(input('Username:'), getpass())
    return auth


def _request_release_data(owner, repo, tag, auth):
    github_api = "https://api.github.com"
    release_endpoint = "/repos/{owner}/{repo}/releases/tags/{tag}".format(
        owner=owner, repo=repo, tag=tag)

    url = github_api + release_endpoint

    logging.critical('Requesting release data: {}'.format(url))

    r = requests.get(url, stream=True, auth=auth)

    if not r.status_code == requests.codes.ok:
        msg = '{}\n{}'.format(url, r.text)
        raise requests.RequestException(msg)
    else:
        return r


def _download_file(url, auth):
    t = tempfile.NamedTemporaryFile()

    logging.critical('Downloading: {}'.format(url))
    r = requests.get(url, stream=True, auth=auth)

    if not r.status_code == requests.codes.ok:
        msg = '{}\n{}'.format(url, r.text)
        raise requests.RequestException(msg)
    else:
        with open(t.name, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
    return t


def _extract_release(path, target):
    """Github release tarballs contain a root directory with the pattern:
    owner-repo-commit and we want them to be referenced by their tag in
    the application. To do this, extract the tarball contents into a temp
    dir, then move all the directory (assuming there should only be one
    directory in the tempdir) to the target path. """
    if os.path.exists(target):
        raise FileExistsError(target)

    # Extract the tarball into a tempdir
    td = tempfile.TemporaryDirectory()
    tf = tarfile.open(path)
    tf.extractall(td.name)

    # Making an assumption here that the only item in the tempdir is
    # the extracted repo.
    items = os.listdir(td.name)
    if len(items) != 1:
        raise RuntimeError('Extracting release: {}, {}'.format(path, items))
    repo_dir = os.path.join(td.name, items[0])
    shutil.move(repo_dir, target)

    td.cleanup()
    return os.path.realpath(target)

# TODO install the latest release
# TODO install all releases
# TODO list all releases

def install_release(tag, owner='tgen', repo='jetstream_pipelines',
                    use_basic_auth=True, credentials=None):
    """ Adds a release to the jetstream pipelines directory. This
    function can actually be used to download any release specified
    with owner/repo. """

    # TODO Handle releases with additional assets
    # Establish auth credentials (needed for private repos)
    if use_basic_auth:
        auth = _basic_auth()
    else:
        auth = credentials

    # Request info about the release
    release_data = _request_release_data(owner, repo, tag, auth)

    # Download the release archive to a tempfile
    archive_url = release_data.json()['tarball_url']
    release_archive_temp = _download_file(archive_url, auth)

    # Extract the release contents to the
    target = os.path.join(RELEASE_DIR, tag)
    return _extract_release(release_archive_temp.name, target)


def list():
    return os.listdir(RELEASE_DIR)


def remove(tag):
    # TODO finish this?
    raise NotImplementedError
