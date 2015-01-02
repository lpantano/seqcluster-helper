import pysam
import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from sqhelper import logger
from sqhelper.group import star_align

import shutil


def quality(data, args):
    _fastqc(data["clean_fastq"])
    return True


def _fastqc(input_file):
    cmd = ("fastqc {input_file}")
    out_file = os.path.abspath(input_file + "_fastqc")
    if not file_exists(out_file):
        do.run(cmd.format(**locals()), "Doing Fastqc")
    return out_file


def _sample_align(data, args):
    work_dir = os.path.join(data["sample_id"], "align")
    safe_makedir(work_dir)
    out_prefix = data["sample_id"] + "_"
    work_dir = os.path.abspath(os.path.join(work_dir, out_prefix))
    out_file = star_align(data, args, data['collapse'], out_prefix)
    return out_file
