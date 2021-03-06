import os
import os.path as op
from bcbio.utils import splitext_plus, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils


def protac(data, args):
    work_dir = data["sample_id"]
    work_dir = op.abspath(utils.safe_makedir(work_dir))
    in_file = data["clean_fastq"]
    with utils.chdir(work_dir):
        collapse = _collapse(in_file)
        clean = _clean(collapse)
        mapper = _mapper(clean, args.reference)
        reallocate = _reallocate(mapper)
        final = _protac(reallocate, args.reference)
    return final


def _collapse(in_file):
    tool = do.find_cmd("TBr2_collapse.pl")
    cmd = "perl {tool} -i {in_file} -o {tx_out}"
    basename = splitext_plus(op.basename(in_file))[0]
    out_file = basename + "_collapse.fastq"
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            do.run(cmd.format(**locals()), "collapse")
    return out_file


def _clean(in_file):
    tool = do.find_cmd("TBr2_duster.pl")
    cmd = "perl {tool} -i {in_file}"
    out_file = in_file + ".no-dust"
    if not utils.file_exists(out_file):
        do.run(cmd.format(**locals()), "duster")
    return out_file


def _mapper(in_file, reference):
    tool = do.find_cmd("sRNAmapper.pl")
    cmd = "perl {tool} -i {in_file} -g {reference} -a best -o {tx_out}"
    out_file = "hits.eland"
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            do.run(cmd.format(**locals()), "mapper")
    return out_file


def _reallocate(in_file):
    tool = do.find_cmd("reallocate.pl")
    cmd = "perl {tool} -i {in_file} 5000 1000 b"
    out_file = in_file + ".weighted-5000-1000-b"
    if not utils.file_exists(out_file):
        do.run(cmd.format(**locals()), "reallocate")
    return out_file


def _protac(in_file, reference):
    tool = do.find_cmd("proTRAC.pl")
    cmd = "perl {tool} -genome  {reference} -map {in_file} -nh -nr -rpm -distr 1-90 -pimax 32"
    out_file = "protac"
    if not utils.file_exists(out_file):
        do.run(cmd.format(**locals()), "protac")
        open(out_file, 'w').close()
    return out_file
