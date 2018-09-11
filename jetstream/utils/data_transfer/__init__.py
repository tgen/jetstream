"""Tools for moving data around TGen.


"""
import requests
from .rsync import rsync, rsync_via_ssh
from .data_movers import Wrangler, RsyncLimitedWrangler, LowestLoadWrangler,\
    check_load, check_rsyncs


def http_download(url, path):
    with open(path, "wb") as fp:
        res = requests.get(url)
        fp.write(res.content)

    return res.ok
