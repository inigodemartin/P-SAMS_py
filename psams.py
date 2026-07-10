#!/usr/bin/env python3

from pathlib import Path
from src.utils import load_config, connect_database, create_opt_subopt, amirna_json, syntasirna_json, parse_list, create_outputs, check_output, get_transcripts, check_accessions, syn_optimal_site_index, load_resume_state
from src.args_parser import parse_args, select_construct
from src.input_parser import convert_fasta_to_string, build_fg_index, build_fg_index_fasta
from src.pipeline import get_tsites, group_tsites, score_sites, serial_jobs
from src.add_species import create_kmer_db
from src.oligo_design import make_amirna_oligos, select_vector, vector_filename_suffix, print_syn_vector_instructions

import os
import shutil
import sys

#############
# CONSTANTS #
#############

PROJECT_ROOT = Path(__file__).resolve().parent
CONFIG_FILE = PROJECT_ROOT / 'psams.conf'
TARGETFINDER = PROJECT_ROOT / 'TargetFinder' / 'targetfinder.pl'
TMP = PROJECT_ROOT / 'tmp'
SEED = 15
MIN_LENGTH = 9



def pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs=1, limit=3, gene_set=None, existing_opt=None, existing_subopt=None, seen_guides=None, start_count=0):

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
    site_scores, TARGETFINDER, mRNA_fa, conn, offtarget, subopt_name, opt_name, output_folder, accession_list, tf_results, len(site_scores), unlimit, fasta, jobs, limit, gene_set,
    existing_opt=existing_opt, existing_subopt=existing_subopt, seen_guides=seen_guides, start_count=start_count)


    ##################################
    # no se como llamar a esta parte #
    ##################################
    opt_count, subopt_count, opt_results, subopt_results = create_opt_subopt(subopt, opt, construct, no_offtarget=not offtarget)

    return opt_count, subopt_count, opt_results, subopt_results


def _warn_if_insufficient(opt_count: int, limit: int, unlimit: bool, gene_set: int = None) -> None:
    """Warn when a run found no optimal sites, or fewer than the requested -l/--limit."""
    label = f" for gene set {gene_set}" if gene_set is not None else ""
    if opt_count == 0:
        print(f"WARNING: no optimal sites were found{label}.")
    elif not unlimit and opt_count < limit:
        print(f"WARNING: only {opt_count} optimal site(s) were found{label} (fewer than the requested limit of {limit}).")


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

    # For syntasiRNA, -V doesn't select a vector here at all: which sites to
    # clone (and the vector) can only be chosen once the optimal sites are
    # known, which happens in the separate clone_vector.py step below.
    vector = None
    if args.vector and construct == "amiRNA":
        vector = select_vector(construct)

    if noofftarget:
        offtarget = False
    else:
        offtarget = True

    # syntasiRNA can define several ':'-separated gene sets, each targeted
    # by its own guide(s); the output folder used to always be auto-named
    # after the first gene set alone, silently ignoring the others in its
    # name. With more than one gene set that's misleading, so a run name is
    # required instead of guessing a folder name for the user. Passing
    # -r/--run_name AND a custom -o/--output_path would be redundant (the
    # last path component of -o already reads as a name), so a custom -o
    # doubles as the run name here: -o runs/gen1_and_gen3 behaves like
    # -o runs -r gen1_and_gen3 used to.
    if construct == "syntasiRNA":
        group_source = fasta if fasta else accessions
        n_gene_sets = len(group_source.split(":"))
        if n_gene_sets > 1 and not args.run_name:
            if output_path != Path("."):
                args.run_name = output_path.name
                output_path = output_path.parent
            else:
                print(
                    f"ERROR: this syntasiRNA run defines {n_gene_sets} gene sets; "
                    "the output folder can't be auto-named after just the first one.\n"
                    "       Pass -r/--run_name, or a named -o/--output_path, to name it explicitly.",
                    file=sys.stderr,
                )
                sys.exit(1)

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
        # Check if transcript inputs are valid (if that ID is in the DB).
        # For syntasiRNA, "accessions" may hold several ':'-separated gene
        # sets (each itself a comma-separated list); flatten across all
        # gene sets for validation, but name the run after the first gene
        # set only (mirrors the -f/--fasta path above).
        transcripts_list = get_transcripts(mRNA_fa)
        accession_groups = [g.split(',') for g in accessions.split(':')]
        check_accessions([acc for group in accession_groups for acc in group], transcripts_list)
        accession_list = accession_groups[0]

    # create outputr folder and check if outputs already exist
    accession_key = args.run_name or '_'.join(accession_list)
    folder_suffix = "_no_offtarget" if noofftarget else ""
    output_folder = output_path / f"{accession_key}_psams_output{folder_suffix}"
    os.makedirs(output_folder, exist_ok=True)

    # amiRNA and syntasiRNA can be run for the same accession(s) into the
    # same output folder: every construct-specific file below is keyed off
    # run_key, not just accession_key, so a syntasiRNA run never collides
    # with (or gets misread as) a previous amiRNA run's cache/results, and
    # vice versa. amiRNA keeps the original, unsuffixed naming so runs from
    # before this fix are still recognized as already-computed.
    run_key = accession_key if construct == "amiRNA" else f"{accession_key}_{construct}"

    # TargetFinder only runs when off-target checking is enabled, so its
    # per-site result folder is pointless (and never written to) in -n mode.
    # tf_results is shared across constructs: amiRNA's site files are never
    # gene-set-prefixed while syntasiRNA's always are (see serial_jobs'
    # gs_prefix), so their filenames never collide even in the same dir.
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
    base_output = cache_dir / f"{run_key}_psams.json"
    visible_output = output_folder / f"{run_key}_psams.json"

    # check if results for given input are already performed
    results_file = output_folder / f"{run_key}_optimal_results.tsv"
    check_output(results_file, accession_list, output_folder, base_output, vector, construct, limit=limit, unlimit=unlimit, want_vector_info=bool(args.vector), run_key=run_key)

    # Resuming: a previous partial run (e.g. capped by a lower -l/--limit)
    # already wrote optimal_results.tsv. Reuse it instead of truncating it
    # via create_outputs — the pipeline below will load_resume_state() from
    # it and skip TargetFinder for every guide already evaluated.
    resume = (not noofftarget) and results_file.exists()

    if noofftarget:
        # No off-target checking means no optimal/suboptimal TSVs: nothing
        # is actually being classified as optimal vs. suboptimal, so those
        # files would only ever contain a header. Only the JSON is written.
        subopt_name = f"{output_folder}/{run_key}_suboptimal_results.tsv"
        opt_name = f"{output_folder}/{run_key}_optimal_results.tsv"
    elif resume:
        subopt_name = f"{output_folder}/{run_key}_suboptimal_results.tsv"
        opt_name = str(results_file)
    else:
        subopt_name, opt_name = create_outputs(run_key, output_folder, syntasirna=(construct == "syntasiRNA"))

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

        if resume:
            existing_opt, existing_subopt, seen_guides, start_count = load_resume_state(tf_results, opt_name, subopt_name)
        else:
            existing_opt, existing_subopt, seen_guides, start_count = [], [], set(), 0

        opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit, existing_opt=existing_opt, existing_subopt=existing_subopt, seen_guides=seen_guides, start_count=start_count)

        _warn_if_insufficient(opt_count, limit, unlimit)

        amirna_json(opt_count, subopt_count, opt_results, subopt_results, base_output, no_offtarget=noofftarget)

        if vector:
            for results in [opt_results, subopt_results]:
                for res in results.values():
                    fwd, rev = make_amirna_oligos(res['guide'], res['star'], vector)
                    res['oligo1'] = fwd
                    res['oligo2'] = rev
            vector_output = output_folder / f"{run_key}_{vector_filename_suffix(vector)}_psams.json"
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
                if resume:
                    existing_opt, existing_subopt, seen_guides, start_count = load_resume_state(tf_results, opt_name, subopt_name, gene_set=g + 1)
                else:
                    existing_opt, existing_subopt, seen_guides, start_count = [], [], set(), 0
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit, gene_set=g + 1, existing_opt=existing_opt, existing_subopt=existing_subopt, seen_guides=seen_guides, start_count=start_count)
                _warn_if_insufficient(opt_count, limit, unlimit, gene_set=g + 1)
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
                if resume:
                    existing_opt, existing_subopt, seen_guides, start_count = load_resume_state(tf_results, opt_name, subopt_name, gene_set=g + 1)
                else:
                    existing_opt, existing_subopt, seen_guides, start_count = [], [], set(), 0
                opt_count, subopt_count, opt_results, subopt_results = pipeline(transcript_dict, foldback, construct, offtarget, unlimit, conn, mRNA_fa, subopt_name, opt_name, output_folder, tf_results, accession_list, fasta, jobs, limit, gene_set=g + 1, existing_opt=existing_opt, existing_subopt=existing_subopt, seen_guides=seen_guides, start_count=start_count)
                _warn_if_insufficient(opt_count, limit, unlimit, gene_set=g + 1)
                groups[g] = {
                "opt": opt_count,
                "sub": subopt_count,
                "opt_r": opt_results,
                "sub_r": subopt_results,
                }

        syntasirna_json(count, groups, base_output, no_offtarget=noofftarget)

        if args.vector:
            # Choosing sites/order/vector always happens as a separate,
            # later step (clone_vector.py) — never inline here — since it's
            # a decision only the user can make and can only be made once
            # these optimal sites are known, which is only now.
            site_index = syn_optimal_site_index(groups, count)
            if not site_index:
                print("No optimal syn-tasiRNA sites were found for any gene set; nothing to clone.")
            else:
                print_syn_vector_instructions(site_index, output_folder, accession_key)

        shutil.copy2(base_output, visible_output)

        print("Finished running P-SAMS successfully!")


if __name__ == "__main__":
    main()