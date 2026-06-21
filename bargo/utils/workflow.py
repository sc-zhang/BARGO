import pathos
import os
from bargo.utils.io import *
from bargo.utils.score import *
from bargo.utils.message import Message


def processing(
    sample,
    bamfile_A,
    bamfile_B,
    gff3_file_A,
    gff3_file_B,
    cds_file_A,
    cds_file_B,
    homo_list_file,
    parent_name_A,
    parent_name_B,
    tau,
    gamma,
    min_mapq,
    out_file,
):
    Message.info("\tProcessing %s..." % sample)
    Message.info("\tLoading GFF3...")
    gff3_A = GFF3(gff3_file_A)
    gff3_B = GFF3(gff3_file_B)

    gff3_A.read_gff3()
    gff3_B.read_gff3()

    Message.info("\tLoading CDS...")
    cds_A = FastaIO(cds_file_A)
    cds_B = FastaIO(cds_file_B)

    cds_A.read_fasta()
    cds_B.read_fasta()

    Message.info("\tProcessing BAM...")
    homolog = Homolog(homo_list_file)
    homolog.read_homolog()

    info = []
    idx = 1

    for geneA, geneB in homolog.homo_pairs:
        reads_A = BamOperator.fetch_gene_reads(bamfile_A, gff3_A.gene_models[geneA])
        reads_B = BamOperator.fetch_gene_reads(bamfile_B, gff3_B.gene_models[geneB])
        L, n_union, n_A, n_B = compute_gene_scores_dual(reads_A, reads_B, min_mapq)
        if L is None:
            gene_class = "Missing"
            conf = float("nan")
            pA = float("nan")
            pB = float("nan")
            L = float("nan")
        else:
            pA, pB = posterior(L)
            conf = confidence(L, n_union)
            gene_class_val = classify(pA, pB, conf, tau, gamma)
            if gene_class_val == 0:
                gene_class = parent_name_A
            elif gene_class_val == 1:
                gene_class = parent_name_B
            else:
                gene_class = "Undetermined"
        if geneA in cds_A.fa_db and geneB in cds_B.fa_db and cds_A.fa_db[geneA] == cds_B.fa_db[geneB]:
            gene_class = "Same"
        info.append([idx, geneA, geneB, n_union, n_A, n_B, L, pA, pB, conf, gene_class])
        idx += 1

    Message.info("\tWriting output file...")
    with open(out_file, "w") as f:
        f.write("#Tau=%f, Gamma=%f\n" % (tau, gamma))
        f.write(
            "#GeneIdx\tRefA\tRefB\tn_union_reads\tn_A_reads\tn_B_reads\tL\tP_A\tP_B\tConfidence\tClass\n"
        )
        for idx, geneA, geneB, n_union, n_A, n_B, L, pA, pB, conf, gene_class in info:
            if gene_class == "Missing":
                f.write(
                    "%d\t%s\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tMissing\n"
                    % (idx, geneA, geneB, n_union, n_A, n_B)
                )
            else:
                f.write(
                    "%d\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%s\n"
                    % (
                        idx,
                        geneA,
                        geneB,
                        n_union,
                        n_A,
                        n_B,
                        L,
                        pA,
                        pB,
                        conf,
                        gene_class,
                    )
                )

    Message.info("\tDone.")


def pipeline(opts):
    parentA = opts.A
    parentB = opts.B
    A_gff3 = opts.A_gff3
    B_gff3 = opts.B_gff3
    A_cds = opts.A_cds
    if not A_cds:
        A_cds = ""
    B_cds = opts.B_cds
    if not B_cds:
        B_cds = ""
    homo_list = opts.list
    tau = opts.tau
    gamma = opts.gamma
    min_mapq = opts.min_mapq
    parent_names = opts.labels.split(",")
    out_dir = opts.output
    threads = opts.threads

    os.makedirs(out_dir, exist_ok=True)

    Message.info("Starting processing...")
    pool = pathos.multiprocessing.Pool(processes=threads)
    res = []
    for fn in os.listdir(parentA):
        if not fn.endswith(".bam"):
            continue
        full_path_A = os.path.join(parentA, fn)
        full_path_B = os.path.join(parentB, fn)
        if not os.path.exists(full_path_B):
            continue
        sample = fn.split(".")[0]
        out_file = os.path.join(out_dir, sample + ".tsv")
        r = pool.apply_async(
            processing,
            (
                sample,
                full_path_A,
                full_path_B,
                A_gff3,
                B_gff3,
                A_cds,
                B_cds,
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
