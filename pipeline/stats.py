import os
import gzip
from counts import _get_position, _get_first_read
from collections import defaultdict


def stats_from_read1(count_file):
    read_gene, counts_gene = _get_first_read(count_file)
    sam_read1 = count_file.replace("assign.dat", "cleaned.sam")
    pos_reads = _get_position(sam_read1, read_gene)
    return pos_reads, read_gene


def stats_from_read2(summary_log_file):
    read2 = defaultdict()
    with open(summary_log_file) as in_handle:
        for line in in_handle:
            if line.find("False") == -1 and line.startswith("corrected"):
                name = line.split(" ")[1]
                new = eval(line.split("--->")[1].split("]")[0]+"]")
                read2[name] = new
    return read2


def quality_read2(polya_log_file, read2):
    with gzip.open(polya_log_file) as in_handle:
        for line in in_handle:
            name, s, e, mod, polya, qpolya, seq, qrest = line.strip().split("\t")
            name = name.replace("@","").split(" ")[0]
            if name in read2:
                read2[name].append(qpolya + qrest)
    return read2


def summarize_stats(count_file, summary_log_file, polya_log_file, out_file):
    gene_pos, read_gene = stats_from_read1(count_file)
    print "read1 done"
    with open(out_file + ".pcr", "w") as out_handle:
        for gene in gene_pos:
            for pos in gene_pos[gene]:
                out_handle.write("%s\n" % " ".join([gene, str(pos),  str(gene_pos[gene][pos]) ]))
    read2 = stats_from_read2(summary_log_file)
    print "read2 done"
    read2 = quality_read2(polya_log_file, read2)
    print "quality done"
    with open(out_file, "w") as out_handle:
        for read in read2:
            if read in read_gene:
                out_handle.write("%s\n" % " ".join([read_gene[read][0], read] + read2[read]))
