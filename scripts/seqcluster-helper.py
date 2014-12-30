from argparse import ArgumentParser
import os
from sqhelper import adapter, align
from sqhelper import cluster


def get_sample(line, sample_map_filename):
    keys = ["sample_id", "fastq", "group"]
    sample_id, r1_filename, group = line.strip().split(",")
    return dict(zip(keys, [sample_id, r1_filename, group]))


def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        return [get_sample(x, sample_file) for x in in_handle]


def write_summary(data):
    with open("summary.csv", 'w') as in_handle:
        for sample in data[0]:
            in_handle.write("%s\n" % ",".join(sample.values()))

if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--minf", type=int,
                        default=2, help="Minimumm counts")
    parser.add_argument("--minl", type=int,
                        default=17, help="Minimumm size")
    parser.add_argument("--db",
                        help="miRBase databases")
    parser.add_argument("--species",
                        default="hsa", help="species")
    parser.add_argument("--aligner", default="STAR", help="STAR or bowtie.")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--aligner-index", help="Path to aligner index.")
    parser.add_argument("--gtf-file", required=False, help="GTF file")
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
    parser.add_argument("--parallel", choices=["local", "ipython"], default="local",
                        help="Run in parallel on a local machine.")
    parser.add_argument("--local", action="store_true",
                        default=False, help="Run parallel locally")

    args = parser.parse_args()

    data = get_samples_to_process(args.sample_map)

    data = cluster.send_job(adapter.remove, data, args, "adapter")
    data = cluster.send_job(align.run_seqcluster, [data], args, "align")

    write_summary(data)
    #cluster.send_job(align.qc, data, args, "qc")

    #data = cluster.send_job(annotate.seqcluster, data, args, "annotate")
