#script to detect modificatino in read2 of tail-seq
import gzip
import re
import numpy
import os
from pandas import Series
from collections import Counter
from utils import open_fastq, file_transaction
from tailseq import logger

w = {}
w['A'] = -8
w['C'] = -8
w['G'] = -8
w['T'] = 2
w['N'] = -5

def getsc(nt):
        if w.has_key(nt):
                return(w[nt])
        else:
                return -10


def poly_A_percentage(seq):
    predicted = {}
    if re.search("TTT", seq):
        for s in range(0, 30):
            #print s
            negative = 0
            scoreT = 0
            if seq[s:].startswith("TT"):
            #string[s]=="T":
                for e in range(s+5, s+60, 2):
                    #print e
                    if negative > 9:
                        break
                    sc = sum([1 for nt in seq[s:e] if nt == "T"])
                    scoreT = 1.0 * sc / (e-s)
                    #print "%s %s %s" % (s, e, scoreT)
                    if scoreT > 0.70:
                        predicted[(e-s+1, s)] = (s, e)
                    else:
                        negative += 1
                if scoreT > 0.95:
                    break
    if len(predicted.keys()) > 0:
        sorted_scores = sorted(predicted.keys(), reverse=True)
        max_ties = [sc[1] for sc in sorted_scores if sc[0] == sorted_scores[0][0]]
        sooner = sorted(max_ties)[0]
        #print "max score %s pos %s " % (sorted_scores[0][0], predicted[(sorted_scores[0][0], sooner)])
        return predicted[(sorted_scores[0][0], sooner)]
    else:
        False


def polyA(string, maxl=30, endl=50):
    smax = -1
    emax = -1
    scmax = -100
    step = 5
    if len(string) < 8:
        step = 1
    if re.search("TTT", string):
        for s in range(0, maxl):
            sc = 0
            converge = 0
            if string[s:].startswith("TT"):
                for e in range(s+step, endl-2, 2):
                    if converge > 10:
                        break
                    if emax >= e - 3 and scmax < 0:
                        converge += 1
                    sc = numpy.sum(map(getsc, string[s:e]))
                    nT = string[s:e].count("T")
                    nO = len(string[s:e])-nT
                    if nO > 0:  # other nt
                        sc = nT - nO + sc
                    if string[e] == "T" and string[e-1] != "T":  # end on NT
                        sc += 2
                        if sc >= scmax+6:
                            scmax = sc
                            smax = s
                            emax = e+1
                            if s > 0:
                                if string[s-1] == "T":
                                    smax = s-1
                                    scmax += 2
                    if string[e-2] == "T" and string[e-1] != "T":  # end on TNN
                        sc += 8
                        if sc >= scmax+6:
                            scmax = sc
                            smax = s
                            emax = e-1
                            if s > 0:
                                if string[s-1] == "T":
                                    smax = s-1
                                    scmax += 2
                    if sc >= scmax + 6 and string[e-1] == "T":  # end on TN
                        scmax = sc
                        smax = s
                        emax = e
                        if s > 0:
                            if string[s-1] == "T":
                                smax = s-1
                                scmax += 2
        if scmax > 0:
            return([smax, emax])
    else:
            return(False)


def _adapter(seq, qual):
    """detect 5nt TAG (GTCAG) and remove from sequence.
    It should be aroung 15-22 position in read
    """
    TAG = "GTCAG"
    tag_pos = seq.find(TAG, 15, 25)
    #print tag_pos
    if tag_pos >= 0:
        tag_pos += 5
        return seq[tag_pos:], qual[tag_pos:]
    else:
        return False


def _test():
        s1 = "CCCCGCATTAAACTTGTCAGAACCAGAGTNATCTTTTTTTTATTTTTTATCTTTTTGATTTATTTTCAGCTCTTCTTTTTCAGTCAAGAATTCTTGCTATTAGGAAAATAATTCCAGATACCATTATAGTAAATATTGCTAAAATGCAAAATACTAATAAAACCTTAGTAAAGTATGAAACTAAAACTAATAGGAAAATTAGAATTGGTGATAATGCTGATAATGAACAATATGAAGTAA"
        s1 = "CAACTGTCGAGGGGGGTCAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTTCGTTGTTTTCTTTTCTATTTATTATTTATTT"
        seq, qual = _adapter(s1, s1)
        print seq
        ns = poly_A_percentage(seq)
        print ns
        print seq[:ns[0]]
        print seq[ns[0]:ns[1]]


def detect(data, args):
    in_file = data['r2_path']
    out_prefix = data['sample_id']
    out_file = out_prefix + "_polyA.dat.gz"
    out_name_false = out_prefix + "_none.dat.gz"
    counts = Counter()
    num_line = 0
    logger.my_logger.info("reading file %s" % in_file)
    logger.my_logger.info("creating files %s %s" % (out_file, out_name_false))
    data['detect'] = out_file
    if os.path.exists(out_file):
        return data
    with file_transaction(out_file) as tx_out_file:
        with open_fastq(in_file) as handle, gzip.open(tx_out_file, 'w') as out, gzip.open(out_name_false, 'w') as out_false:
            for line in handle:
                #print line
                num_line += 1
                if num_line % 1000000 == 0:
                    logger.my_logger.info("read %s lines:" % num_line)
                if line.startswith("@HISEQ"):
                    #print line
                    name = line.strip()
                    seq = handle.next().strip()
                    handle.next().strip()
                    qual = handle.next().strip()
                    find = _adapter(seq, qual)
                    #print "%s %s" % (seq, find)
                    if find:
                        seq, qual = find
                        ns = poly_A_percentage(seq)
                        #ns = polyA(seq)
                        if ns:
                            if ns[1]-ns[0] >= 6:
                                #print "positions are" + str(ns[0]) + ".." + str(ns[1])
                                mod = seq[:ns[0]]
                                seq_polyA = seq[ns[0]:ns[1]]
                                seq_gene = seq[ns[1]:]
                                qual_polyA = qual[ns[0]:ns[1]]
                                qual_gene = qual[ns[1]:]
                                #print "%s\t%s\t%s\t%s\t%s\t%s\n" % (name,mod,sf,qf)
                                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, ns[0], ns[1], mod, seq_polyA, qual_polyA, seq_gene, qual_gene))
                                counts['polyA'] += 1
                                if len(mod) > 0:
                                    counts['mod'] += 1
                            else:
                                counts['shortA'] += 1
                                out_false.write("%s\t%s\t%s\t%s\n" % ("shortA", name, seq, qual))
                        else:
                            counts['noA'] += 1
                            out_false.write("%s\t%s\t%s\t%s\n" % ("None", name, seq, qual))
                    else:
                        out_false.write("%s\t%s\t%s\t%s\n" % ("No_tag", name, seq, qual))
                        counts['notag'] += 1
        with file_transaction(out_prefix + ".stat") as tx_stat_file:
            df = Series(counts)
            df.to_csv(tx_stat_file, sep="\t")
        logger.my_logger.info("%s" % counts)
    return data


def tune(mod, polya):
    """correct mod"""
    if sum([1 for nt in mod if nt == "T"]) == len(mod):
        return ["", mod+polya]
    if len(polya) < 5:
        return False
    seq = mod + polya
    start_limit = len(mod) + 15 if len(polya) > 25 else len(mod) + len(polya) - 4
    end_limit = 50 + len(mod) if len(polya) > 54 else len(seq)
    find = polyA(seq, start_limit, end_limit)
    if find:
        if len(polya) > 4:
            mod = seq[:find[0]]
            polya = seq[find[0]:find[1]]
            return [mod, polya]
        return False
    return False
