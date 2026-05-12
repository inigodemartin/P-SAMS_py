#!/usr/bin/env python3

from pathlib import Path
import argparse
from src.utils import index_genome, add_species_to_config
from src.add_species import create_kmer_db, add_species
from src.create_sqlite_db import create_sqlite_db
from src.phytozome2psams import phytozome2psams
import os 


PROJECT_ROOT = Path(__file__).resolve().parent

def parse_args():
    parser = argparse.ArgumentParser(
        description="Adds a new species to the kmer database (Python version)."
    )
    parser.add_argument("-f", "--fasta", required=True, help="FASTA transcript file")
    parser.add_argument("-s", "--species", required=True, help="Species code")
    parser.add_argument("-k", "--ksize", required=True, type=int, help="Kmer size")
    parser.add_argument("-v", "--version", help="Transcriptome version")
    parser.add_argument("-d", "--descriptions", help="Annotation description file")
    parser.add_argument("-l", "--lowmem", action="store_true", help="Enable low memory mode")
    return parser.parse_args()



def main():
    args = parse_args()
    fasta = Path(args.fasta)

    # index genome
    indexed_fasta = fasta.parent / f"{fasta.name}.fai"
    if indexed_fasta.exists() and indexed_fasta.stat().st_size > 0:
        pass    
    else:
        print(f"Indexing {fasta}") 
        index_genome(fasta)

    # Change phytozome format
    
    formatted_fasta, formatted_annot = phytozome2psams(fasta, args.species, args.version, args.descriptions, fasta.parent)
    db_file = PROJECT_ROOT / 'db' / f"{args.species}_{args.version}_prueba.db"

    if args.lowmem:
        # create kmers and add directly to the db, using low RAM but much more time
        print(f"Creating {args.species} DB in low memory mode...")
        create_kmer_db(db_file, formatted_fasta, args.ksize, formatted_annot)
    else:
        # create kmer file
        print("Creating kmer file")
        kmer_file = add_species(formatted_fasta, args.ksize)

        # load kamer file and descriptions to the database
        print("Creating database")
        
        create_sqlite_db(db_file, kmer_file, formatted_annot)

    # Load database info in config file
    add_species_to_config(PROJECT_ROOT / "psams.conf", args.species, db_file, fasta)



if __name__ == "__main__":
    main()