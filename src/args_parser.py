import argparse

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

    return args
