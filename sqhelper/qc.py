import pysam
import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from sqhelper import logger
from sqhelper.group import star_align

import shutil

def quality(data, args):
    out_prefix = data["sample_id"] + "_"
    star_align(data, args, data['collapse'], out_prefix)
    return True
