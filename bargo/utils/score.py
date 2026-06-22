import numpy as np
from collections import defaultdict


def alignment_score(read, min_mapq=20):
    mapq = read.mapping_quality
    if mapq < min_mapq:
        return None
    as_score = read.get_tag("AS") if read.has_tag("AS") else 0
    return as_score


def compute_gene_scores_dual(reads_a, reads_b, min_mapq=20):
    read_dict = defaultdict(lambda: [0.0, 0.0])
    read_set_a = set()
    read_set_b = set()
    for read in reads_a:
        qname = read.query_name
        s = alignment_score(read, min_mapq)
        if s is None:
            continue
        read_dict[qname][0] = s
        read_set_a.add(qname)

    for read in reads_b:
        qname = read.query_name
        s = alignment_score(read, min_mapq)
        if s is None:
            continue
        read_dict[qname][1] = s
        read_set_b.add(qname)
    deltas = []
    for sa, sb in read_dict.values():
        deltas.append(sa - sb)
    if len(deltas) == 0:
        return None, 0, 0, 0
    L = np.mean(deltas)
    return L, len(deltas), len(read_set_a), len(read_set_b)


def posterior(L):
    if L >= 0:
        z = np.exp(-L)
        p_a = 1.0 / (1.0 + z)
    else:
        z = np.exp(L)
        p_a = z / (1.0 + z)
    p_b = 1.0 - p_a
    return p_a, p_b


def confidence(L, n_reads):
    return abs(L) * (n_reads / (n_reads + 5))


# 0: class A, 1: class B, 2: undetermined
def classify(p_a, p_b, conf, tau=0.7, gamma=0.3):
    if conf < gamma:
        return 2
    if p_a > tau:
        return 0
    if p_b > tau:
        return 1
    return 2
