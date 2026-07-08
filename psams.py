#!/usr/bin/env python3

from pathlib import Path
from src.utils import load_config, connect_database, create_opt_subopt, amirna_json, syntasirna_json, parse_list, create_outputs, check_output, get_transcripts, check_accessions
from src.args_parser import parse_args, select_construct
from src.input_parser import convert_fasta_to_string, build_fg_index, build_fg_index_fasta
from src.pipeline import get_tsites, group_tsites, score_sites, serial_jobs
from src.add_species import create_kmer_db
from src.oligo_design import make_amirna_oligos, make_syntasirna_oligos, select_vector, prompt_target_site, vector_filename_suffix

import os
import shutil

#############
# CONSTANTS #
#############

PROJECT_ROOT = Path(__file__).resolve().parent
CONFIG_FILE = PROJECT_ROOT / 'psams.conf'
TARGETFINDER = PROJECT_ROOT / 'TargetFinder' / 'targetfinder.pl'
TMP = PROJECT_ROOT / 'tmp'
SEED = 15
MIN_LENGTH = 9



def pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs=1, limit=3):

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
    site_scores, TARGETFINDER, mRNA_fa, conn, offtarget, subopt_name, opt_name, output_folder, accession_list, tf_results, len(site_scores), unlimit, fasta, jobs, limit)


    ##################################
    # no se como llamar a esta parte #
    ##################################
    opt_count, subopt_count, opt_results, subopt_results = create_opt_subopt(subopt, opt, construct, no_offtarget=not offtarget)

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
    construct = args.construct or select_construct()
    unlimit = args.unlimit
    limit = args.limit
    jobs = args.jobs
    noofftarget = args.noofftarget
    output_path = Path(args.output_path)
    target_site = args.target_site

    vector = None
    if args.vector:
        vector = select_vector(construct)
        if vector == "pMDC32B-B/c" and not target_site:
            target_site = prompt_target_site()

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
    accession_key = '_'.join(accession_list)
    folder_suffix = "_no_offtarget" if noofftarget else ""
    output_folder = output_path / f"{accession_key}_psams_output{folder_suffix}"
    os.makedirs(output_folder, exist_ok=True)

    # TargetFinder only runs when off-target checking is enabled, so its
    # per-site result folder is pointless (and never written to) in -n mode.
    tf_results = output_folder / "tf_results"
    if not noofftarget:
        os.makedirs(tf_results, exist_ok=True)

    # base_output is the canonical, vector-agnostic result cache: it is
    # always (re)written on a full run and reused by check_output below to
    # skip the expensive pipeline when only the cloning vector changes. It
    # lives in a hidden .cache/ subfolder so it never shows up next to the
    # vector-specific result and gets mistaken for a second, valid result.
    cache_dir = output_folder / ".cache"
    os.makedirs(cache_dir, exist_ok=True)
    base_output = cache_dir / f"{accession_key}_psams.json"
    visible_output = output_folder / f"{accession_key}_psams.json"

    # check if results for given input are already performed
    results_file = output_folder / f"{accession_key}_optimal_results.tsv"
    check_output(results_file, accession_list, output_folder, base_output, vector, target_site, construct)

    if noofftarget:
        # No off-target checking means no optimal/suboptimal TSVs: nothing
        # is actually being classified as optimal vs. suboptimal, so those
        # files would only ever contain a header. Only the JSON is written.
        subopt_name = f"{output_folder}/{accession_key}_suboptimal_results.tsv"
        opt_name = f"{output_folder}/{accession_key}_optimal_results.tsv"
    else:
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

        opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit)

        amirna_json(opt_count, subopt_count, opt_results, subopt_results, base_output, no_offtarget=noofftarget)

        if vector:
            for results in [opt_results, subopt_results]:
                for res in results.values():
                    fwd, rev = make_amirna_oligos(res['guide'], res['star'], vector)
                    res['oligo1'] = fwd
                    res['oligo2'] = rev
            vector_output = output_folder / f"{accession_key}_{vector_filename_suffix(vector)}_psams.json"
            amirna_json(opt_count, subopt_count, opt_results, subopt_results, vector_output, vector=vector, no_offtarget=noofftarget)
        else:
            # No vector selected: the cache is the only result, surface it.
            shutil.copy2(base_output, visible_output)

        print("Finished running P-SAMS successfully!")

    elif construct == "syntasiRNA":
        
        groups = {}
        count = 0

        if fasta:
            fasta_groups = fasta.split(":")
            count = len(fasta_groups)

            for g in range(count):
                fasta_str = convert_fasta_to_string(fasta_groups[g])
                transcript_dict = build_fg_index_fasta(fasta_str)
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit)
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
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit)
                groups[g] = {
                "opt": opt_count,
                "sub": subopt_count,
                "opt_r": opt_results,
                "sub_r": subopt_results,
                }

        syntasirna_json(count, groups, base_output, no_offtarget=noofftarget)

        if vector:
            syn_guides = []
            for g in range(count):
                group = groups[g]
                if group["opt"] > 0:
                    syn_guides.append(group["opt_r"][1]["guide"])
                elif group["sub"] > 0:
                    syn_guides.append(group["sub_r"][1]["guide"])
            cloning_oligos = None
            if syn_guides:
                fwd, rev = make_syntasirna_oligos(syn_guides, vector, target_site)
                cloning_oligos = {"forward": fwd, "reverse": rev}
            vector_output = output_folder / f"{accession_key}_{vector_filename_suffix(vector)}_psams.json"
            syntasirna_json(count, groups, vector_output, vector=vector, cloning_oligos=cloning_oligos, no_offtarget=noofftarget)
        else:
            # No vector selected: the cache is the only result, surface it.
            shutil.copy2(base_output, visible_output)

        print("Finished running P-SMAS successfully!")


if __name__ == "__main__":
    main()