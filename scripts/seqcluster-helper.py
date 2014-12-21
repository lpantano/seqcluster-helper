from argparse import ArgumentParser
import os
from tailseq import align
from tailseq import cluster
from tailseq import detect
from tailseq import counts


def get_sample(line, sample_map_filename):
    keys = ["sample_id", "r1_path", "r2_path"]
    sample_id, r1_filename, r2_filename = line.strip().split(",")
    return dict(zip(keys, [sample_id, r1_filename, r2_filename]))


def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        return [get_sample(x, sample_file) for x in in_handle]


def get_r2_prepped_outfile(sample, alignment_dir):
    return os.path.join(alignment_dir,
                        ".".join([sample["sample_id"], sample["subsample_id"]]))


def get_star_prefix(fastq_file):
    base, _ = os.path.splitext(fastq_file)
    return base + "_base"


def get_cleaned_outfile(align_file):
    base, ext = os.path.splitext(align_file)
    return base + "_cleaned" + ext


def get_prefix(in_file):
    return in_file.split("_base")[0]


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--multimappers", action="store_true",
                        default=False, help="Keep multimappers")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--aligner-index", help="Path to aligner index.")
    parser.add_argument("--alignment-dir", help="Output directory")
    parser.add_argument("--gtf-file", required=True, help="GTF file")
    parser.add_argument("--num-jobs", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--scheduler", default=None, help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge"])
    parser.add_argument("--resources", default=None, help="Extra scheduler resource flags.")
    parser.add_argument("--queue", default=None, help="Queue to submit jobs to.")
    parser.add_argument("--parallel", choices = ["local", "ipython"], default="local",
                        help="Run in parallel on a local machine.")
    parser.add_argument("--local", action="store_true",
                        default=False, help="Run parallel locally")
    parser.add_argument("--rmdup", action="store_true",
                        default=False, help="remove duplicates")

    args = parser.parse_args()

    data = get_samples_to_process(args.sample_map)

    data = cluster.send_job(align.align_read1, data, args, "align_r1")

    data = cluster.send_job(align.align_read2, data, args, "align_r2")

    cluster.send_job(align.qc, data, args, "qc")

    data = cluster.send_job(detect.detect, data, args, "polyA")

    data = cluster.send_job(counts.counts, data, args, "counts")
