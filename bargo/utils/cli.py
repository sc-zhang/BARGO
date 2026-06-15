import argparse
from bargo.__version__ import __version__
from bargo.utils.workflow import pipeline


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-A", help="Input bam directory of parent A", required=True)
    group.add_argument("-B", help="Input bam directory of parent B", required=True)
    group.add_argument("--A_gff3", help="Input gff3 file of parent A", required=True)
    group.add_argument("--B_gff3", help="Input gff3 file of parent B", required=True)
    group.add_argument(
        "-l", "--list", help="Input homologous gene pair list", required=True
    )
    group.add_argument(
        "--labels",
        help='Name of parents, comma split, default="A,B"',
        type=str,
        default="A,B",
    )
    group.add_argument(
        "--tau",
        help="Posterior decision threshold, default=0.7",
        type=float,
        default=0.7,
    )
    group.add_argument(
        "--gamma",
        help="Confidence filtering threshold, default=0.3",
        type=float,
        default=0.3,
    )
    group.add_argument(
        "--min_mapq", help="Minimum mapping quality, default=20", type=float, default=20
    )
    group.add_argument(
        "-o",
        "--output",
        help="Output directory of tsv files of classification",
        required=True,
    )
    group.add_argument(
        "-t",
        "--threads",
        help="Number of threads for parallel processing, default=1",
        default=1,
        type=int,
    )
    group.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return group.parse_args()


def main():
    opts = get_opts()
    pipeline(opts)
