import pathos
import os
from bargo.utils.io import *
from bargo.utils.score import *
from bargo.utils.message import Message


def processing(
    sample,
    bamfile_a,
    bamfile_b,
    gff3_file_a,
    gff3_file_b,
    identical_cds_pair_file,
    homo_list_file,
    parent_name_a,
    parent_name_b,
    tau,
    gamma,
    min_mapq,
    out_file,
):
    Message.info("\tProcessing %s..." % sample)
    Message.info("\tLoading GFF3...")
    gff3_a = GFF3(gff3_file_a)
    gff3_b = GFF3(gff3_file_b)

    gff3_a.read_gff3()
    gff3_b.read_gff3()

    Message.info("\tLoading identical CDS pairs...")
    identical_cds = GenePair(identical_cds_pair_file)
    identical_cds.read_pairs()
    identical_cds_set = set()
    for geneA, geneB in identical_cds.pairs:
        identical_cds_set.add(tuple([geneA, geneB]))

    Message.info("\tProcessing BAM...")
    homolog = GenePair(homo_list_file)
    homolog.read_pairs()

    info = []
    idx = 1

    for geneA, geneB in homolog.pairs:
        reads_a = BamOperator.fetch_gene_reads(bamfile_a, gff3_a.gene_models[geneA])
        reads_b = BamOperator.fetch_gene_reads(bamfile_b, gff3_b.gene_models[geneB])
        L, n_union, n_a, n_b = compute_gene_scores_dual(reads_a, reads_b, min_mapq)
        if L is None:
            gene_class = "Missing"
            conf = float("nan")
            p_a = float("nan")
            p_b = float("nan")
            L = float("nan")
        else:
            p_a, p_b = posterior(L)
            conf = confidence(L, n_union)
            gene_class_val = classify(p_a, p_b, conf, tau, gamma)
            if gene_class_val == 0:
                gene_class = parent_name_a
            elif gene_class_val == 1:
                gene_class = parent_name_b
            else:
                gene_class = "Undetermined"
        if tuple([geneA, geneB]) in identical_cds_set:
            gene_class = "Identical"
        info.append(
            [idx, geneA, geneB, n_union, n_a, n_b, L, p_a, p_b, conf, gene_class]
        )
        idx += 1

    Message.info("\tWriting output file...")
    with open(out_file, "w") as f:
        f.write("#Tau=%f, Gamma=%f\n" % (tau, gamma))
        f.write(
            "#GeneIdx\tRefA\tRefB\tn_union_reads\tn_A_reads\tn_B_reads\tL\tP_A\tP_B\tConfidence\tClass\n"
        )
        for idx, geneA, geneB, n_union, n_a, n_b, L, p_a, p_b, conf, gene_class in info:
            if gene_class == "Missing":
                f.write(
                    "%d\t%s\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tMissing\n"
                    % (idx, geneA, geneB, n_union, n_a, n_b)
                )
            else:
                f.write(
                    "%d\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%s\n"
                    % (
                        idx,
                        geneA,
                        geneB,
                        n_union,
                        n_a,
                        n_b,
                        L,
                        p_a,
                        p_b,
                        conf,
                        gene_class,
                    )
                )

    Message.info("\tDone.")


def pipeline(opts):
    parent_a = opts.A
    parent_b = opts.B
    a_gff3_file = opts.A_gff3
    b_gff3_file = opts.B_gff3
    a_cds_file = opts.A_cds
    if not a_cds_file:
        a_cds_file = ""
    b_cds_file = opts.B_cds
    if not b_cds_file:
        b_cds_file = ""
    homo_list = opts.list
    tau = opts.tau
    gamma = opts.gamma
    min_mapq = opts.min_mapq
    parent_names = opts.labels.split(",")
    out_dir = opts.output
    threads = opts.threads

    os.makedirs(out_dir, exist_ok=True)

    Message.info("Starting processing...")

    Message.info("\tPreprocessing...")
    identical_cds_pair_file = ""

    if a_cds_file and b_cds_file:
        cds_a = FastaIO(a_cds_file)
        cds_b = FastaIO(b_cds_file)

        cds_a.read_fasta()
        cds_b.read_fasta()

        homolog = GenePair(homo_list)
        homolog.read_pairs()

        identical_cds_pair_file = os.path.join(out_dir, "identical_cds_pairs.tsv")
        with open(identical_cds_pair_file, "w") as f:
            for geneA, geneB in homolog.pairs:
                if (
                    geneA in cds_a.fa_db
                    and geneB in cds_b.fa_db
                    and cds_a.fa_db[geneA] == cds_b.fa_db[geneB]
                ):
                    f.write("%s\t%s\n" % (geneA, geneB))

    result_dir = os.path.join(out_dir, "result")
    os.makedirs(result_dir, exist_ok=True)

    pool = pathos.multiprocessing.Pool(processes=threads)
    res = []
    for fn in os.listdir(parent_a):
        if not fn.endswith(".bam"):
            continue
        full_path_a = os.path.join(parent_a, fn)
        full_path_b = os.path.join(parent_b, fn)
        if not os.path.exists(full_path_b):
            continue
        sample = fn.split(".")[0]
        out_file = os.path.join(result_dir, sample + ".tsv")
        r = pool.apply_async(
            processing,
            (
                sample,
                full_path_a,
                full_path_b,
                a_gff3_file,
                b_gff3_file,
                identical_cds_pair_file,
                homo_list,
                parent_names[0],
                parent_names[1],
                tau,
                gamma,
                min_mapq,
                out_file,
            ),
        )
        res.append(r)
    pool.close()
    pool.join()

    for r in res:
        try:
            r.get()
        except Exception as e:
            Message.error(repr(e))

    Message.info("All done.")
