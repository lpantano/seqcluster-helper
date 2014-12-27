import os
from collections import Counter
from bcbio.utils import splitext_plus, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils
import shutil

ADAPTER = "AGATCGGAAGAGCAC"


def remove(data, args):
    adapter = ADAPTER
    work_dir = data["sample_id"]
    work_dir = os.path.abspath(utils.safe_makedir(work_dir))
    in_file = data["fastq"]
    out_dir = os.path.join(work_dir, "adapter")
    utils.safe_makedir(out_dir)
    basename = splitext_plus(os.path.basename(in_file))[0]
    out_file = os.path.join(out_dir, basename + "_clean.fastq")
    out_noadapter_file = os.path.join(out_dir, basename + "_fragments.fastq")
    out_short_file = os.path.join(out_dir, basename + "_short.fastq")
    cmd = _cmd_cutadapt()
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "remove adapter")
    data["clean_fastq"] = out_file
    data['collapse'] = _collapse(out_file)
    data['size_stats'] = _summary(data['collapse'])
    out_dir = utils.safe_makedir(os.path.join(work_dir, 'miraligner'))
    out_file = os.path.join(out_dir, data["sample_id"])
    data['miraligner'] = _miraligner(data["clean_fastq"], out_file,args. species, args.db)
    return data


def _cmd_cutadapt():
    cmd = "cutadapt --adapter={adapter} --minimum-length=8 --untrimmed-output={out_noadapter_file} -o {tx_out_file} -m 17 --overlap=8 {in_file} --too-short-output {out_short}"
    return cmd


def _collapse(in_file):
    cmd = "seqcluster collapse -f {in_file} -o {out_dir}"
    basename = splitext_plus(os.path.basename(in_file))[0]
    out_file = splitext_plus(in_file)[0] + "_trimmed.fastq"
    if not utils.file_exists(out_file):
        with tx_tmpdir() as out_dir:
            tx_out_file = os.path.join(out_dir, basename + "_trimmed.fastq")
            do.run(cmd.format(**locals()), "collapse")
            shutil.move(tx_out_file, out_file)
    return out_file


def _summary(in_file):
    data = Counter()
    out_file = in_file + "_size_stats"
    with open(in_file) as in_handle:
        for line in in_handle:
            counts = int(line.strip().split("_x")[1])
            line = in_handle.next()
            l = len(line.strip())
            in_handle.next()
            in_handle.next()
            data[l] += counts
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for l, c in data.iteritems():
                out_handle.write("%s %s\n" % (l, c))
    return out_file


def _miraligner(fastq_file, out_file, species, db_folder):
    cmd = ("miraligner -Xms750m -Xmx8g -sub 1 -trim 3 -add 3 -s {species} -i {fastq_file} -db {db_folder}  -o {tx_out_file}")
    if not file_exists(out_file + ".mirna"):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "Do miRNA annotation")
            shutil.move(tx_out_file + ".mirna", out_file + ".mirna")
    return out_file
