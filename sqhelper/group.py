import pysam
import os
from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.log import logger
import shutil

MAX_EDIT_DISTANCE = 2
MAX_BEST = 1000

logger.LOG_NAME = "seqcluster-helper"

def run_seqcluster(data, args):
    out_dir = "seqcluster"
    minl, minf = args.minl, args.minf
    out_dir = os.path.abspath(safe_makedir(out_dir))
    config_file = os.path.join(out_dir, "prepare.conf")
    prepare_dir = os.path.join(out_dir, "prepare")
    prepare_dir = _prepare(data, config_file, prepare_dir, minl, minf)
    fastq_file = os.path.join(prepare_dir, "seqs.fastq")
    bam_file = _align(data, fastq_file, args)
    #quality = _coverage(bam_file, prepare_dir)
    cluster_dir = os.path.join(out_dir, "cluster")
    cluster_dir = _cluster(bam_file, prepare_dir, cluster_dir, args.reference, args.gtf_file)
    data = _update(data, prepare_dir, bam_file, cluster_dir)
    return data


def _prepare(data, config_file, out_dir, minl, minf):
    with file_transaction(config_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for sample in data:
                fasta = sample['collapse']
                name = sample['sample_id']
                out_handle.write("%s\t%s\n" % (fasta, name))
    min_samples = min(0, int(round(len(data) / 10)))
    cmd = ("seqcluster prepare -c {config_file} -u 40 -l {minl} -e {minf} -o {tx_out_dir} --min-shared {min_samples}")
    if not file_exists(out_dir):
        with tx_tmpdir() as work_dir:
            tx_out_dir = os.path.join(work_dir, "prepare")
            do.run(cmd.format(**locals()), "seqcluster prepare")
            shutil.move(tx_out_dir, out_dir)
    return out_dir


def _align(data, fastq_file, args):
    work_dir = os.path.join("align")
    work_dir = os.path.abspath(safe_makedir(work_dir))
    out_prefix = os.path.join(work_dir, "seqs_")
    bam_file = star_align(data, args, fastq_file, out_prefix, 1000)
    return bam_file


def _cluster(bam_file, prepare_dir, out_dir, reference, annotation_file="None"):
    opts = ""
    if annotation_file:
        opts = "-g %s" % annotation_file
    cmd = ("seqcluster cluster -m {ma_file} -ref {reference} -a {bam_file} -o {out_dir} {opts}")
    ma_file = os.path.join(prepare_dir, "seqs.ma")
    if not file_exists(os.path.join(out_dir, "counts.tsv")):
        # out_dir = os.path.join(out_dir, "cluster")
        do.run(cmd.format(**locals()), "seqcluster cluster")
    return out_dir


def star_align(data, args, fastq_path, out_prefix, max_hits=10, opts=""):
    cores = args.cores_per_job
    reference_prefix = args.aligner_index
    out_file = out_prefix + "Aligned.sortedByCoord.out.bam"
    clean_file = out_prefix + "clean.bam"
    if not os.path.exists(out_file):
        cmd = ("STAR --genomeDir {reference_prefix} --readFilesIn {fastq_path} --outSAMtype BAM SortedByCoordinate "
           "--runThreadN {cores} --outFileNamePrefix {out_prefix} "
           "--outFilterMultimapNmax {max_hits} "
           "--outSAMattributes NH HI NM {opts} "
           "--alignIntronMax 1 ").format(**locals())
        do.run(cmd, "Alignment")
    clean_align(out_file, clean_file)
    return clean_file


def _coverage(bam_file, prepare_dir):
    ma_file = os.path.join(prepare_dir, "seqs.ma")
    list_seq = parse_ma_file(bam_file)
    with pysam.Samfile(ma_file, "rb") as bam:
        for a in bam.fetch():
            if not a.is_unmapped:
                if a.qname in list_seq:
                    for sample in list_seq:
                        sam_files[sample] += a


def _update(data,  prepare_dir, bam_file, cluster_dir):
    for s in data:
        s["seqs_ma"] = os.path.join(prepare_dir, "seqs.ma")
        s["seqs_bam"] = bam_file
        s["clusters"] = os.path.join(cluster_dir, "counts.tsv")
    return data


def qc(data, args):
    """fastqc for the sam file"""
    sam_file = data['align']
    out_dir = os.path.basename(sam_file) + "_fastq"
    cmd = "fastqc {sam_file} -f sam -o {out_dir}".format(**locals())
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        do.run(cmd)
    else:
        logger.info("%s has already been QC, skipping." % (sam_file))
    return data


def clean_align(align_file, out_file):
    if file_exists(out_file):
        logger.info("%s has already been cleaned, skipping." % (align_file))
        return out_file

    count_total_reads = 0
    count_assigned_reads = 0
    count_assigned_aligned_reads = 0
    with pysam.Samfile(align_file, "rb") as in_handle, file_transaction(out_file) as tx_out_file:
        out_handle = pysam.Samfile(tx_out_file, "wb", template=in_handle)
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
    nm = get_tag(read, "NM")
    if nm and (nm > MAX_EDIT_DISTANCE):
        return True
    splice = [True for cigar in read.cigar if cigar[0] == 3]
    if any(splice):
        return True
#    x0 = get_tag(read, "X0")
#    if x0 and (x0 > MAX_BEST):
#        return True
    return False


def get_tag(read, tag):
    matched_tag = [x for x in read.tags if tag == x[0]]
    return matched_tag[0][1] if matched_tag else None
