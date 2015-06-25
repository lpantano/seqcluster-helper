import os.path as op
from argparse import ArgumentParser
from sqhelper import sample, group, qc
from sqhelper import cluster
from sqhelper import pirna
from sqhelper.report import create_rmd
from sqhelper import do


def get_sample(line, sample_map_filename):
    keys = ["fastq", "sample_id", "group"]
    if len(line.split(',')) < 3:
        raise ValueError("This line hasn't at least 3 elements (name,path_to_fasta,group): %s" % line)
    cols = line.strip().split(",")
    r1_filename, sample_id, group = cols[:3]
    assert op.exists(r1_filename), "File doesn't exists: %s" % r1_filename
    return dict(zip(keys, [r1_filename, sample_id, group]))


def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        in_handle.next()
        return [get_sample(x, sample_file) for x in in_handle]


def write_summary(data):
    header = data[0][0].keys()
    with open("summary.csv", 'w') as in_handle:
        in_handle.write("%s\n" % ",".join(header))
        for s in data[0]:
            in_handle.write("%s\n" % ",".join(s.values()))
    return "summary.csv"


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
    parser.add_argument("--adapter", required=True, help="Adapter.")
    parser.add_argument("--aligner-index", help="Path to aligner index.")
    parser.add_argument("--reference", help="Path to genome fasta.")

    parser.add_argument("--gtf-file", required=False, help="GTF file")
    parser.add_argument("--protac", action='store_true', required=False, help="Run protac pipeline.")
    parser.add_argument("--config", required=False, help="yaml file with custom resources.")
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
    parser.add_argument("--report", action="store_true",
                        default=False, help="Run Rmarkdown to create html")

    args = parser.parse_args()

    data = get_samples_to_process(args.sample_map)

    data = cluster.send_job(sample.remove, data, args, "sample")
    data = cluster.send_job(qc.quality, data, args, "qc")

    data = cluster.send_job(group.run_seqcluster, [data], args, "group")

    summary_file = write_summary(data)
    rel_summary_file, report_file = create_rmd(summary_file)

    if args.report:
        do.run('R -e "library(rmarkdown);library(knitrBootstrap);render(\'%s\')"' % report_file, "Create html")

    print "your summary files and report are: %s and %s\n" % (rel_summary_file, report_file)
    print "you can open %s with R/Rstudio and create the report\n" % report_file

    if args.protac:
        cluster.send_job(pirna.protac, [data], args, "pirna")
