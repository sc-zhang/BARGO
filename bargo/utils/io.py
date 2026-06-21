import pysam


class FastaIO:
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        self.fa_db = {}

    def read_fasta(self):
        if self.fasta_file != "":
            with open(self.fasta_file, "r") as f:
                for line in f:
                    if line[0] == ">":
                        gid = line.strip().split()[0][1:]
                        self.fa_db[gid] = []
                    else:
                        self.fa_db[gid].append(line.strip())

            for gid in self.fa_db:
                self.fa_db[gid] = "".join(self.fa_db[gid])


class GeneModel:
    def __init__(self, chrom, cds_regions: list):
        self.chrom = chrom
        self.cds_regions = cds_regions


class GFF3:
    def __init__(self, gff3_file: str):
        self.gff3_file = gff3_file
        self.gene_models = {}

    def read_gff3(self):
        with open(self.gff3_file, "r") as f:
            for line in f:
                if line.strip() == "" or line[0] == "#":
                    continue
                data = line.strip().split("\t")
                chrom = data[0]
                feature_type = data[2]
                if feature_type == "gene":
                    gene_id = ""
                    gname = ""
                    for info in data[8].split(";"):
                        if info.startswith("ID"):
                            gene_id = info.split("=")[1]
                        if info.startswith("Name"):
                            gname = info.split("=")[1]
                    if gname != "":
                        gene_id = gname
                    if gene_id not in self.gene_models:
                        self.gene_models[gene_id] = GeneModel(chrom, [])
                elif feature_type == "CDS":
                    start = int(data[3])
                    end = int(data[4])
                    self.gene_models[gene_id].cds_regions.append([start, end])


class Homolog:
    def __init__(self, homolog_file: str):
        self.homolog_file = homolog_file
        self.homo_pairs = []

    def read_homolog(self):
        with open(self.homolog_file, "r") as f:
            for line in f:
                data = line.strip().split("\t")
                self.homo_pairs.append([data[0], data[1]])


class BamOperator:
    @staticmethod
    def fetch_gene_reads(bam_file, gene: GeneModel):
        reads = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            seen = set()
            for start, end in gene.cds_regions:
                for read in bam.fetch(gene.chrom, start, end):
                    if read.is_unmapped:
                        continue
                    qname = read.query_name
                    if qname in seen:
                        continue
                    seen.add(qname)
                    reads.append(read)
        return reads
