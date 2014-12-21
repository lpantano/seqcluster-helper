from utils import file_transaction
import os
import pysam
import gzip
from collections import defaultdict, Counter
from tailseq import do
from tailseq import logger
from tailseq.detect import tune


def counts(data, args):
    logger.my_logger.info("Counting sample %s:" % data['sample_id'])
    in_file = data['clean']
    prefix = data['sample_id'] + "_"
    gtf_file = args.gtf_file
    cores = args.cores_per_job
    out_file = prefix + "counts"
    data['counts'] = _cmd_counts(in_file, out_file, gtf_file, cores)
    data['assign'] = _assign_gene(data['counts'], prefix)
    _summarize(data['detect'], data['align_r2'], data['assign'], prefix + "summary.dat")
    return data


def _cmd_counts(in_file, out_file, gtf_file, cores):
    if not os.path.exists(out_file):
        cmd = "featureCounts -R -T {cores} --primary -a {gtf_file} -o {out_file} {in_file}"
        do.run(cmd.format(**locals()))
    return in_file + ".featureCounts"


def _assign_gene(in_file, prefix):
    """read featureCounts output and assign each read a gene"""
    out_file = prefix + "assign.dat"
    if not os.path.exists(out_file):
        with open(in_file) as handle, file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out:
                for line in handle:
                    cols = line.strip().split("\t")
                    if cols[1] == "Assigned":
                        out.write("%s\t%s\n" % (cols[0], cols[2]))
    return out_file


def _summarize(in_file, align_r2, count_file, out_file):
    log_file = out_file + ".log"
    logger.my_logger.info("summarize results")
    read_gene, counts_gene = _get_first_read(count_file)
    logger.my_logger.info("load read 1 done")
    read_gene = _get_second_read(align_r2, read_gene)
    logger.my_logger.info("load read 2 done")
    stats = defaultdict(Counter)
    if not os.path.exists(out_file):
        with gzip.open(in_file) as handle_polya:
            log_handle = open(log_file, 'w')
            for line in handle_polya:
                cols = line.strip().split("\t")
                read = cols[0].split(" ")[0].replace("@", "")
                if read in read_gene:
                    log_handle.write("found %s %s %s ---> %s\n" % (read, cols[3], cols[4], read_gene[read]))
                    if read_gene[read][1]:
                        if len(cols[3] + cols[4] + cols[6]) > 135:
                            continue
                        find = tune(cols[3], cols[4])
                        if find:
                            log_handle.write("corrected %s %s %s --->%s %s\n" % (read, cols[3], cols[4], find, read_gene[read]))
                            #print "is polya"
                            gene = read_gene[read][0]
                            polya_size = _get_bin(len(find[1]))
                            stats[gene][polya_size] += 1
                            if find[0] != "":
                                stats[gene][(polya_size, find[0])] += 1
                        else:
                            log_handle.write("removed %s %s %s ---> %s %s\n" % (read, cols[3], cols[4], find, read_gene[read]))
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out:
                for gene in counts_gene:
                    out.write("%s total counts %s 0\n" % (gene, counts_gene[gene]))
                    if gene in stats:
                        for mod, c in stats[gene].iteritems():
                            if isinstance(mod, tuple):
                                u_times = _get_u_times(mod[1])
                                out.write("%s %s %s %s %s\n" % (gene, mod[0], mod[1], c, u_times))
                            else:
                                out.write("%s polyA %s %s 0\n" % (gene, mod, c))


def _get_u_times(m):
    us = sum([1 for nt in m if nt == "A"])
    ratio = 1.0 * us / len(m)
    return ratio


def _get_bin(size):
    if size < 15:
        return "<15"
    elif size < 25:
        return "<25"
    else:
        return ">25"


def _get_second_read(sam_file, read1_assing):
    num_lines = 0
    with pysam.Samfile(sam_file, 'r') as sam:
        for read in sam.fetch():
            num_lines += 1
            if num_lines % 10000000 == 0:
                logger.my_logger.info(num_lines)
            if not read.is_unmapped and read.qname in read1_assing and not read.is_secondary:
                if not read.is_read2:
                    continue
                nm = int([t[1] for t in read.tags if t[0] == "NM"][0])
                name = read.qname
                cigar = read.cigar[::-1] if read.is_reverse else read.cigar
                read1_assing[name][1] = _is_polyA(cigar, nm, read.is_proper_pair)
                # if name == "HISEQ:272:HAVVLADXX:1:1101:3750:2216":
                #    print [read.is_reverse, cigar, read.is_read2]
                #    raise
                #print "%s %s %s %s %s %s %s" % (read.seq, name, nm, cigar, read.is_reverse, read.is_proper_pair, _is_polyA(cigar, nm, read.is_proper_pair))
    return read1_assing


def _is_polyA(cigar, nm, pair):
    """guess if all nt after 20 adapter are mapped"""
    if cigar[0][1] == 20 and _get_middle(cigar) > 90 and pair:
        return False
    if cigar[0][1] - 20 >= 6 or _get_middle(cigar) < 50:
        return [cigar, nm]
    return False


def _get_middle(cigar):
    l = 0
    for c in cigar:
        if c[0] != 4:
            l += c[1]
        # else:
        #    break
    return l


def _get_first_read(in_file):
    genes = defaultdict(list)
    stats = Counter()
    with open(in_file) as counts:
        for line in counts:
            cols = line.strip().split("\t")
            genes[cols[0]] = [cols[1], "Added"]
            stats[cols[1]] += 1
    return genes, stats


def _get_position(sam_file, genes):
    """get all positions from SAM file and sync with _get_first_read"""
    pos = defaultdict(Counter)
    with pysam.Samfile(sam_file, 'r') as sam:
        for read in sam.fetch():
            if not read.is_unmapped and read.qname in genes and not read.is_secondary:
                name = read.qname
                pos[genes[name][0]][int(read.pos)] += 1
    return pos
