import os
from bcbio.utils import splitext_plus
from bcbio.distributed import do

ADAPTER = "AATGATACGGCGA"


def remove(data, args):
    work_dir = data["sample_id"]
    in_file = data["fastq"]
    out_dir = os.path.join(work_dir, "adapter")
    out_file = os.path.join(out_dir, splitext_plus(os.path.basename(in_file))[0] + "_clean.fastq")
    out_noadapter_file = os.path.join(out_dir, splitext_plus(os.path.basename(in_file))[0] + "_fragments.fastq")
    cmd = _cmd_cutadapt()
    do.run(cmd.format(**locals()))
    data["clean_fastq"] = out_file
    return data


def _cmd_cutadapt():
    cmd = "cutadapt --adapter={adapter} --minimum-length=8 --untrimmed-output={out_noadapter_file} -o {out_file} -m 17 --overlap=8 {in_file}"
    return cmd
