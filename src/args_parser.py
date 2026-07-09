import argparse


def select_construct() -> str:
    """Interactively prompt the user to choose a construct type."""
    options = ["amiRNA", "syntasiRNA"]

    print("\nSelect construct type:")
    for i, name in enumerate(options, start=1):
        print(f"  {i}. {name}")

    while True:
        choice = input(f"Select construct [1-{len(options)}]: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(options):
            return options[int(choice) - 1]
        print(f"Invalid selection '{choice}'. Please enter a number between 1 and {len(options)}.")


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
    parser.add_argument("-c", "--construct",choices=["amiRNA", "syntasiRNA"],help="Construct type. If omitted, you'll be prompted to choose interactively.")
    parser.add_argument("-n", "--noofftarget",action="store_true",help="Run without predicting off-target transcripts.")
    parser.add_argument("-u", "--unlimit",action="store_true",help="Unlimited results (slow).")
    parser.add_argument("-l", "--limit",type=int,default=3,help="Maximum number of optimal results to find before stopping. Default = 3. Ignored if -u/--unlimit is set.")
    parser.add_argument("-j", "--jobs",type=int,default=1,help="Number of TargetFinder jobs to run in parallel. Default = 1 (serial).")
    parser.add_argument("-V", "--vector",
        action="store_true",
        help="For amiRNA, prompts to interactively pick a cloning vector and adds vector-specific "
             "cloning oligos to the output. For syntasiRNA, prints the optimal sites found and the "
             "exact clone_vector.py command to run afterwards to generate the cloning oligos — "
             "choosing sites/vector/order always happens in that separate step, once the sites are "
             "known (see clone_vector.py --help).")
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

    if args.limit < 1:
        parser.error("-l/--limit must be at least 1.")

    return args
