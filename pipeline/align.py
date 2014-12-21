import re
from utils import file_transaction, file_exists
import pysam
import os
from tailseq import do
from tailseq import logger

accepted_barcode_pattern = re.compile(r"[ACGT]+[ACG]$")
polyA_tail_pattern = re.compile(r"A{20,}$")
MAX_EDIT_DISTANCE = 15
MAX_BEST = 10


def align_read1(data, args):
    fastq_path = data['r1_path']
    out_prefix = data['sample_id'] + "_"
    data["align"], data["clean"] = star_align(data, args, fastq_path, out_prefix)
    if args.rmdup:
        out_file = data['sample_id'] + "_rmdup.bam"
        data["clean"] = rmdup(data["clean"], out_file)
    return data


def align_read2(data, args):
    fastq_path = data['r1_path'] + " " + data['r2_path']
    out_prefix = data['sample_id'] + "_both_"
    opts = "--clip3pNbases 0 25 --clip5pNbases 0 20 "
    data["align_r2"], data["clean_r2"] = star_align(data, args, fastq_path, out_prefix, opts)
    return data


def star_align(data, args, fastq_path, out_prefix, opts=""):
    cores = args.cores_per_job
    reference_prefix = args.aligner_index
    max_best = MAX_BEST
    out_file = out_prefix + "Aligned.out.sam"
    if not os.path.exists(out_file):
        cmd = ("STAR --genomeDir {reference_prefix} --readFilesIn {fastq_path} --readFilesCommand zcat "
           "--runThreadN {cores} --outFileNamePrefix {out_prefix} "
           "--outFilterMultimapNmax {max_best} "
           "--outSAMattributes NH HI NM MD AS {opts} "
           "").format(**locals())
        do.run(cmd)
    clean_file = clean_align(out_file, out_prefix + "cleaned.sam")
    return out_file, clean_file


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


def rmdup(align_file, out_file):
    cmd = ("samtools view -Sbh {align_file} | samtools sort -o -n - {tmp} | bammarkduplicates rmdup=1   O={tx_out_file}")
    tmp = align_file + "_tmp"
    if not os.path.exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()))
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
