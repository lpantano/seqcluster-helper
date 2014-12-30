import pysam
import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from sqhelper import logger
import shutil

MAX_EDIT_DISTANCE = 15
MAX_BEST = 1000


def run_seqcluster(data, args):
    out_dir = "seqcluster"
    out_dir = os.path.abspath(safe_makedir(out_dir))
    config_file = os.path.join(out_dir, "prepare.conf")
    prepare_dir = os.path.join(out_dir, "prepare")
    prepare_dir = _prepare(data, config_file, prepare_dir)
    fastq_file = os.path.join(prepare_dir, "seqs.fastq")
    bam_file = _align(data, fastq_file, args)
    #quality = _coverage(bam_file, prepare_dir)
    cluster_dir = os.path.join(out_dir, "cluster")
    cluster_dir = _cluster(bam_file, prepare_dir, cluster_dir, args.gtf_file)
    return data


def _prepare(data, config_file, out_dir):
    with file_transaction(config_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for sample in data:
                fasta = sample['collapse']
                name = sample['sample_id']
                out_handle.write("%s\t%s\n" % (fasta, name))
    cmd = ("seqcluster prepare -c {config_file} -u 40 -l 16 -e 5 -o {tx_out_dir}")
    if not file_exists(out_dir):
        with tx_tmpdir() as work_dir:
            tx_out_dir = os.path.join(work_dir, "prepare")
            do.run(cmd.format(**locals()), "seqcluster prepare")
            shutil.move(tx_out_dir, out_dir)
    return out_dir


def _align(data, fastq_file, args):
    work_dir = os.path.join("align")
    work_dir = os.path.abspath(safe_makedir(work_dir))
    out_prefix = os.path.join(work_dir, "seqs")
    bam_file = star_align(data, args, fastq_file, out_prefix)
    return bam_file


def _cluster(bam_file, prepare_dir, out_dir, annotation_file="None"):
    if annotation_file:
        opts = "-b %s" % annotation_file
    cmd = ("seqcluster cluster -m {ma_file} -a {bam_file} -o {tx_out_dir} {opts} -d ")
    ma_file = os.path.join(prepare_dir, "seqs.ma")
    if not file_exists(out_dir):
        with tx_tmpdir() as work_dir:
            tx_out_dir = os.path.join(work_dir, "cluster")
            do.run(cmd.format(**locals()), "seqcluster cluster")
            shutil.move(tx_out_dir, out_dir)
    return out_dir


def star_align(data, args, fastq_path, out_prefix, opts=""):
    cores = args.cores_per_job
    reference_prefix = args.aligner_index
    max_best = MAX_BEST
    out_file = out_prefix + "Aligned.sortedByCoord.out.bam"
    if not os.path.exists(out_file):
        cmd = ("STAR --genomeDir {reference_prefix} --readFilesIn {fastq_path} --outSAMtype BAM SortedByCoordinate "
           "--runThreadN {cores} --outFileNamePrefix {out_prefix} "
           "--outFilterMultimapNmax {max_best} "
           "--outSAMattributes NH HI NM {opts} "
           "").format(**locals())
        do.run(cmd, "Alignment")
    return out_file


def _coverage(bam_file, prepare_dir):
    ma_file = os.path.join(prepare_dir, "seqs.ma")
    list_seq = parse_ma_file(bam_file)
    with pysam.Samfile(ma_file, "rb") as bam:
        for a in bam.fetch():
            if not a.is_unmapped:
                if a.qname in list_seq:
                    for sample in list_seq:
                        sam_files[sample] += a


def qc(data, args):
    """fastqc for the sam file"""
    sam_file = data['align']
    out_dir = os.path.basename(sam_file) + "_fastq"
    cmd = "fastqc {sam_file} -f sam -o {out_dir}".format(**locals())
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        do.run(cmd)
    else:
        logger.my_logger.info("%s has already been QC, skipping." % (sam_file))
    return data


def clean_align(align_file, out_file):
    if file_exists(out_file):
        logger.my_logger.info("%s has already been cleaned, skipping." % (align_file))
        return out_file

    count_total_reads = 0
    count_assigned_reads = 0
    count_assigned_aligned_reads = 0
    with pysam.Samfile(align_file, "r") as in_handle, file_transaction(out_file) as tx_out_file:
        out_handle = pysam.Samfile(tx_out_file, "wh", template=in_handle)
        for read in in_handle:
            count_total_reads += 1
            count_assigned_reads += 1
            if poorly_mapped_read(read):
                continue
            count_assigned_aligned_reads += 1
            out_handle.write(read)
        out_handle.close

    return out_file


def poorly_mapped_read(read):
    if read.is_unmapped:
        return True
#    nm = get_tag(read, "NM")
#    if nm and (nm > MAX_EDIT_DISTANCE):
#        return True
    x0 = get_tag(read, "X0")
    if x0 and (x0 > MAX_BEST):
        return True
    return False


def get_tag(read, tag):
    matched_tag = [x for x in read.tags if tag == x[0]]
    return matched_tag[0][1] if matched_tag else None
