import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.bam import fastq
from sqhelper import logger
from sqhelper.group import star_align

# import shutil


def quality(data, args):
    work_dir = os.path.join(data["sample_id"], "qc")
    data["qc"] = _fastqc(data["clean_fastq"], work_dir)
    # data["align_alone"] = _sample_align(data, args)
    return data


def _fastqc(input_file, out_dir):
    data = {'config': {'algorithm': {}}}
    if not file_exists(out_dir):
        print input_file
        dw_file, _ = fastq.downsample(input_file, None, data, int(1e7))
        cmd = ("fastqc {dw_file} --extract -o {out_dir}")
        out_dir = os.path.abspath(out_dir)
        safe_makedir(out_dir)
        do.run(cmd.format(**locals()), "Doing Fastqc %s" % input_file)
        logger.my_logger.debug(cmd.format(**locals()))
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
