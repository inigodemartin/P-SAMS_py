import argparse
from src.oligo_design import AMIRNA_VECTORS, SYNTASIRNA_VECTORS

_ALL_VECTORS = list(AMIRNA_VECTORS) + list(SYNTASIRNA_VECTORS)

_AMIRNA_VECTOR_HELP = "\n  ".join(AMIRNA_VECTORS)
_SYNTASIRNA_VECTOR_HELP = "\n  ".join(SYNTASIRNA_VECTORS)


def parse_args():

    """
    Parse command-line arguments for the P-SAMS-like tool.
    """
    parser = argparse.ArgumentParser(description="Plant Small RNA Maker Suite (P-SAMS). Artificial miRNA and syn-tasiRNA designer tool.",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--fasta",help="FASTA-formatted sequence. Not used if -a is set.")
    parser.add_argument("-a", "--accessions",help="Gene accession(s). Comma-separated list. Not used if -f is set.")
    parser.add_argument("-s", "--species",help="Species. Required if -a is set.")
    parser.add_argument("-o", "--output_path",default=".",help="Output folder path.")
    parser.add_argument("-t", "--foldback",default="eudicot",choices=["eudicot", "monocot"],help="Foldback type [eudicot, monocot]. Default = eudicot.")
    parser.add_argument("-c", "--construct",default="amiRNA",choices=["amiRNA", "syntasiRNA"],help="Construct type. Default = amiRNA.")
    parser.add_argument("-n", "--noofftarget",action="store_true",help="Run without predicting off-target transcripts.")
    parser.add_argument("-u", "--unlimit",action="store_true",help="Unlimited results (slow).")
    parser.add_argument("-j", "--jobs",type=int,default=1,help="Number of TargetFinder jobs to run in parallel. Default = 1 (serial).")
    parser.add_argument("-V", "--vector",
        choices=_ALL_VECTORS,
        metavar="VECTOR",
        help=(
            "Cloning vector for oligo design. Generates vector-specific cloning oligos in the output.\n"
            "amiRNA vectors:\n"
            f"  {_AMIRNA_VECTOR_HELP}\n"
            "syn-tasiRNA vectors:\n"
            f"  {_SYNTASIRNA_VECTOR_HELP}"
        ))
    parser.add_argument("-T", "--target-site",
        dest="target_site",
        metavar="SEQ",
        help="22-nt miRNA target site sequence. Required when --vector pMDC32B-B/c is used with --construct syntasiRNA.")
    parser.add_argument("-p", "--phytozome_fasta", help="Genome assembly fasta")
    parser.add_argument("-d", "--descriptions", help="Annotation description file")
    parser.add_argument("-v", "--version", help="Genome version")

    args = parser.parse_args()

    if not args.accessions and not args.fasta:
        parser.error("An input sequence or a gene accession ID were not provided!")

    if args.accessions and not args.species:
        parser.error("A species name was not provided!")

    if args.fasta and not args.species and not args.noofftarget:
        parser.error("Off-target prediction with -f requires -s/--species (or use -n/--noofftarget to disable it).")

    if args.jobs < 1:
        parser.error("-j/--jobs must be at least 1.")

    if args.vector:
        if args.construct == "amiRNA" and args.vector not in AMIRNA_VECTORS:
            parser.error(
                f"Vector '{args.vector}' is a syn-tasiRNA vector and cannot be used with --construct amiRNA.\n"
                f"Available amiRNA vectors: {', '.join(AMIRNA_VECTORS)}"
            )
        if args.construct == "syntasiRNA" and args.vector not in SYNTASIRNA_VECTORS:
            parser.error(
                f"Vector '{args.vector}' is an amiRNA vector and cannot be used with --construct syntasiRNA.\n"
                f"Available syn-tasiRNA vectors: {', '.join(SYNTASIRNA_VECTORS)}"
            )
        if args.vector == "pMDC32B-B/c" and not args.target_site:
            parser.error("--target-site is required when using --vector pMDC32B-B/c.")

    if args.target_site:
        ts = args.target_site.strip().upper().replace("U", "T")
        invalid = set(ts) - set("ACGT")
        if invalid or len(ts) != 22:
            parser.error("--target-site must be a 22-nt DNA sequence (only A, C, G, T bases allowed).")
        args.target_site = ts

    return args
