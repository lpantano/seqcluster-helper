import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
# from sqhelper import logger
from sqhelper.group import star_align

# import shutil


def quality(data, args):
    work_dir = os.path.join(data["sample_id"], "qc")
    data["qc"] = _fastqc(data["clean_fastq"], work_dir)
    # data["align_alone"] = _sample_align(data, args)
    return data


def _fastqc(input_file, out_dir):
    cmd = ("fastqc {input_file} --extract -o {out_dir}")
    out_dir = os.path.abspath(out_dir)
    safe_makedir(out_dir)
    if not file_exists(out_dir):
        do.run(cmd.format(**locals()), "Doing Fastqc")
    return out_dir


def _sample_align(data, args):
    work_dir = os.path.join(data["sample_id"], "align")
    safe_makedir(work_dir)
    out_prefix = data["sample_id"] + "_"
    out_prefix = os.path.abspath(os.path.join(work_dir, out_prefix))
    out_file = star_align(data, args, data['collapse'], out_prefix, 10)
    return out_file


def _coverage(bam_file, out_file):
    cmd = ("bamtools coverage -in {bam_file} -out {tx_out_file}")
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "Coverage of %s" % bam_file)
    return out_file
