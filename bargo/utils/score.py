import numpy as np
from collections import defaultdict


def alignment_score(read, min_mapq=20):
    mapq = read.mapping_quality
    if mapq < min_mapq:
        return None
    as_score = read.get_tag("AS") if read.has_tag("AS") else 0
    return as_score


def compute_gene_scores_dual(reads_A, reads_B, min_mapq=20):
    read_dict = defaultdict(lambda: [0.0, 0.0])
    for read in reads_A:
        qname = read.query_name
        s = alignment_score(read, min_mapq)
        if s is None:
            continue
        read_dict[qname][0] = s

    for read in reads_B:
        qname = read.query_name
        s = alignment_score(read, min_mapq)
        if s is None:
            continue
        read_dict[qname][1] = s
    deltas = []
    for sa, sb in read_dict.values():
        deltas.append(sa - sb)
    if len(deltas) == 0:
        return 0.0, 0
    L = np.mean(deltas)
    return L, len(deltas)


def posterior(L):
    if L >= 0:
        z = np.exp(-L)
        pA = 1.0 / (1.0 + z)
    else:
        z = np.exp(L)
        pA = z / (1.0 + z)
    pB = 1.0 - pA
    return pA, pB


def confidence(L, n_reads):
    return abs(L) * (n_reads / (n_reads + 5))


# 0: class A, 1: class B, 2: undetermined
def classify(pA, pB, conf, tau=0.7, gamma=0.3):
    if conf < gamma:
        return 2
    if pA > tau:
        return 0
    if pB > tau:
        return 1
    return 2
