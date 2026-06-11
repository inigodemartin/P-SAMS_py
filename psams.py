#!/usr/bin/env python3

from pathlib import Path
from src.utils import load_config, connect_database, create_opt_subopt, amirna_json, syntasirna_json, parse_list, create_outputs, check_output, get_transcripts, check_accessions
from src.args_parser import parse_args
from src.input_parser import convert_fasta_to_string, build_fg_index, build_fg_index_fasta
from src.pipeline import get_tsites, group_tsites, score_sites, serial_jobs
from src.add_species import create_kmer_db

import os

#############
# CONSTANTS #
#############

PROJECT_ROOT = Path(__file__).resolve().parent
CONFIG_FILE = PROJECT_ROOT / 'psams.conf'
TARGETFINDER = PROJECT_ROOT / 'TargetFinder' / 'targetfinder.pl'
TMP = PROJECT_ROOT / 'tmp'
SEED = 15
MIN_LENGTH = 9



def pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta):

    #################################
    # Identify putative targe sites #
    #################################

    t_sites = get_tsites(transcript_dict, SEED, offtarget, conn=conn) 
    gsites = group_tsites(t_sites, SEED)

    #################################
    # Score putative targe sites    #
    #################################
    target_count = len(transcript_dict)
    site_scores = score_sites(gsites, target_count, SEED, foldback, debug=False)
    
    ######################################
    # Run TargetFinder on top candidates #
    ######################################
    opt, subopt = serial_jobs(
    target_count, construct, transcript_dict,
    site_scores, TARGETFINDER, mRNA_fa, conn, offtarget, subopt_name, opt_name, output_folder, accession_list, tf_results, len(site_scores), unlimit, fasta)


    ##################################
    # no se como llamar a esta parte #
    ##################################
    opt_count, subopt_count, opt_results, subopt_results = create_opt_subopt(subopt, opt, construct)

    return opt_count, subopt_count, opt_results, subopt_results

def main():
    
    ################################
    # Parse command-line arguments #
    ################################

    args = parse_args()
    
    fasta = args.fasta
    accessions = args.accessions
    species = args.species
    foldback = args.foldback
    construct = args.construct
    unlimit = args.unlimit
    noofftarget = args.noofftarget
    output_path = Path(args.output_path)

    if noofftarget:
        offtarget = False
    else:
        offtarget = True

    #######################################
    # Load config and connect to database #
    #######################################

    config = load_config(CONFIG_FILE)
    conn, mRNA_fa = connect_database(config, species)

    if fasta:
        # FASTA input: derive the naming list from the sequence IDs in the
        # (first) FASTA file, no database lookup needed.
        first_fasta = fasta.split(":")[0]
        accession_list = list(build_fg_index_fasta(convert_fasta_to_string(first_fasta)).keys())
    else:
        # Check if transcript inputs are valid (if that ID is in the DB)
        transcripts_list = get_transcripts(mRNA_fa)
        accession_list = accessions.split(',')
        check_accessions(accession_list, transcripts_list)

    # create outputr folder and check if outputs already exist
    output_folder = output_path / f"{'_'.join(accession_list)}_psams_output"
    tf_results = output_folder / "tf_results"
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(tf_results, exist_ok=True)

    # check if results for given input are already performed
    results_file = output_folder / f"{'_'.join(accession_list)}_optimal_results.tsv"
    check_output(results_file, accession_list, output_folder)

    subopt_name, opt_name = create_outputs(accession_list, output_folder)

    ##################################################################
    # Run the pipeline distinguishing between amiRNA and syntasiRNA  #
    ##################################################################
    
    if construct == "amiRNA":
        if fasta:
            fasta_str = convert_fasta_to_string(fasta)
            transcript_dict = build_fg_index_fasta(fasta_str)

        elif accessions:
            accession_list = accessions.split(',')
            transcript_dict = build_fg_index(accession_list, conn, species, mRNA_fa)

        opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta)
        

        final_output = output_folder / f"{'_'.join(accession_list)}_psams.json"
        amirna_json(opt_count, subopt_count, opt_results, subopt_results, final_output)
        print("Finished running P-SMAS successfully!")

    elif construct == "syntasiRNA":
        
        groups = {}
        count = 0

        if fasta:
            fasta_groups = fasta.split(":")
            count = len(fasta_groups)

            for g in range(count):
                fasta_str = convert_fasta_to_string(fasta_groups[g])
                transcript_dict = build_fg_index_fasta(fasta_str)
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta)
                groups[g] = {
                "opt": opt_count,
                "sub": subopt_count,
                "opt_r": opt_results,
                "sub_r": subopt_results,
                }

        else:

            accession_groups = accessions.split(":")
            count = len(accession_groups)
            for g in range(count):
                accession_list = accession_groups[g].split(',')
                transcript_dict = build_fg_index(accession_list, conn, species, mRNA_fa)
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta)
                groups[g] = {
                "opt": opt_count,
                "sub": subopt_count,
                "opt_r": opt_results,
                "sub_r": subopt_results,
                }

        final_output = output_folder / f"{'_'.join(accession_list)}_psams.json"        
        syntasirna_json(count, groups, final_output)
        
        print("Finished running P-SMAS successfully!")


if __name__ == "__main__":
    main()