import pysam
import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from sqhelper import logger
import shutil

def quality(data, args):
    return True
