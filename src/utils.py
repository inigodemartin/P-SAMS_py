import configparser
from pprint import pprint
import json
import sys
import shutil
import subprocess
import time
import threading
import sqlite3
from pathlib import Path

from src.oligo_design import make_amirna_oligos, make_syntasirna_oligos, vector_filename_suffix, prompt_syn_order

########################################################
#                   CONFIG LOADER                      #
########################################################

def load_config(config_file):
    """
    Load configuration file
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    return config


########################################################
#               DATABASE CONNECTION                    #
########################################################

def connect_database(config, species):
    """
    Connect to SQLite database using config values.
    """

    if not species:
        return None, None

    mRNA_fa = config[species]['mRNA']
    db_path = config[species]['sql']

    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    return connection, mRNA_fa


def design_guide_RNA(site: dict) -> str:
    """
    Design guide RNA based on target site sequences.

    Parameters
    ----------
    site : dict
        Contains 'seqs' (semicolon-separated sequences)

    Returns
    -------
    str
        Guide RNA sequence (reversed)
    """
    guide = ""

    mm = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'AC': 'A', 'AG': 'A', 'AT': 'C',
        'CG': 'A', 'CT': 'C',
        'ACG': 'A', 'ACT': 'C',
        'GT': 'G',
        'AGT': 'G',
        'CGT': 'T',
        'ACGT': 'A'
    }

    bp = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'AC': 'T', 'AG': 'T', 'AT': 'T',
        'CG': 'C', 'CT': 'G',
        'ACG': 'T', 'ACT': 'G',
        'GT': 'C',
        'AGT': 'T',
        'CGT': 'A',
        'ACGT': 'T'
    }

    sites = site['seqs'].split(';')
    # Iterate over 21 positions
    for i in range(21):
        nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # Count nucleotides at position i
        for seq in sites:
            nts[seq[i]] += 1

        # Build diversity string
        nt_str = ''.join(nt for nt in ['A', 'C', 'G', 'T'] if nts[nt] > 0)

        # Apply rules
        if i == 0:
            guide += mm[nt_str] # AQUI SE METE UN MISMATCH INTENCIONAL PARA QUE EL amiRNA QUE SE DISEÑA NO SE PEGUE DEMASIADO (no se exactamente por que)
        elif i == 2:
            guide += 'C' # PARECE SER QUE INTERESA AQUI UNA UNION G-C PQ ACTUA COMO ANCLA INICIAL EN EL SEED, PARA ASEGURAR QUE LA UNIÓN DEL amiRNA SEA ESTABLE DESDE EL INICIO

        elif i == 20:
            guide += 'T' # AQUI SI O SI TIENE QUE SER UNA T PQ LOS MIRNA SIEMPRE LLEVAN UNA T AHI, LO NECESITAN PARA ENGANCHARSE AL RISC
                         # TODOS LOS miRNAS NATURALES TIENEN UN NUCLEÓTIDO U/T EN ESA POSICIÓN PARA INTRODUCIRSE EN EL RISC

        else:
            guide += bp[nt_str] #PARA LAS DEMAS POSICIONES SE ELIGE LA COMPLEMENTARIA DEL TARGET

    # Reverse sequence
    return guide[::-1]

def get_transcripts(fasta):
    transcripts = []
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                id = line.strip().split()[0][1:]
                transcripts.append(id)
    return transcripts

def check_accessions(accession_list, transcript_list):
    for accession in accession_list:
        if accession not in transcript_list:
            sys.exit(f"Error: {accession} not found in database.")


def index_genome(fasta):
    """
    Index a FASTA file using samtools faidx.

    Parameters
    ----------
    fasta : str
        Path to FASTA file.
    """

    cmd = ["samtools", "faidx", fasta]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"samtools faidx failed:\n{result.stderr}"
        )

    return f"{fasta}.fai"



def oligo_designer(guide: str, fb_type: str):
    """
    Generate cloning oligonucleotide sequences.

    Parameters
    ----------
    guide : str
        Guide RNA sequence
    fb_type : str
        'eudicot' or 'monocot'

    Returns
    -------
    tuple
        (realstar, oligo1, oligo2)
    """

    # Reverse complement
    rev = guide[::-1]
    rev = rev.translate(str.maketrans("ACGTacgt", "TGCAtgca"))

    temp = list(rev)
    c = temp[10]
    g = temp[2]
    n = temp[20]

    # Complement of single base 
    c = c.translate(str.maketrans("AGCT", "CTAG"))
    if fb_type == 'eudicot':

        star1 = rev[:10] + c + rev[11:21]

        oligo1 = guide + 'ATGATGATCACATTCGTTATCTATTTTTT' + star1

        oligo2 = oligo1[::-1]
        oligo2 = oligo2.translate(str.maketrans("ACTGacgt", "TGACtgca"))

        realstar = star1[2:22] + 'CA'

        string = 'AGTAGAGAAGAATCTGTA' + oligo1 + 'CATTGGCTCTTCTTACT'
        bsa1 = 'TGTA'
        bsa2 = 'AATG'

    elif fb_type == 'monocot':

        star1 = rev[:10] + c + rev[11:20] + 'C'

        oligo1 = guide + 'ATGATGATCACATTCGTTATCTATTTTTT' + star1

        oligo2 = oligo1[::-1]
        oligo2 = oligo2.translate(str.maketrans("ATGCatgc", "TACGtacg"))

        realstar = star1[2:22] + 'CA'

        string = 'GGTATGGAACAATCCTTG' + oligo1 + 'CATGGTTTGTTCTTACC'
        bsa1 = 'CTTG'
        bsa2 = 'CATG'

    else:
        raise ValueError(f"Foldback type {fb_type} not supported.")

    return realstar, bsa1 + oligo1, bsa2 + oligo2




def _status_loop(start_time, accession_list, optimal_ref, suboptimal_ref, stop_event, potential_target_n, offtarget=True):
    """
    Live status line updated every second.
    Runs in background thread.
    """

    while not stop_event.is_set():

        elapsed = int(time.time() - start_time)

        hours = elapsed // 3600
        minutes = (elapsed % 3600) // 60
        seconds = elapsed % 60

        if offtarget:
            line = (
                f"Running P-SAMS for {','.join(accession_list)} | "
                f"{hours:02d}:{minutes:02d}:{seconds:02d} | "
                f"Optimal sites found: {optimal_ref[0]} | "
                f"Suboptimal sites found: {suboptimal_ref[0]} | "
                f"Potential target sites: {potential_target_n} | "
            )
        else:
            line = f"Running P-SAMS for {','.join(accession_list)} in no-offtarget mode"

        # Truncate to terminal width and clear the rest of the line, so a
        # long line never wraps (which would turn each update into a new
        # line on screen instead of overwriting the current one).
        width = shutil.get_terminal_size((80, 20)).columns
        sys.stdout.write("\r" + line[:width - 1] + "\033[K")
        sys.stdout.flush()

        time.sleep(1)


def create_outputs(accessions, output_folder, syntasirna=False):
    subopt_name = f"{output_folder}/{'_'.join(accessions)}_suboptimal_results.tsv"
    opt_name = f"{output_folder}/{'_'.join(accessions)}_optimal_results.tsv"

    with open(subopt_name, 'w') as subopt, open(opt_name, 'w') as opt:
        if syntasirna:
            # syntasiRNA runs may cover several gene sets writing into this
            # same file (one call per gene set); Gene_set identifies which
            # gene set a row belongs to, so results can be addressed as
            # "geneset.site" (e.g. "1.1", "3.2") when choosing what to clone.
            subopt_header = 'Gene_set\tSite_index\tOfftarget_N\tOfftarget_list\tGuide\tStar\tOligo1\tOligo2\tIsoforms\tIsoform_seqs\n'
            opt_header = 'Gene_set\tSite_index\tGuide\tStar\tOligo1\tOligo2\tIsoforms\tIsoform_seqs\n'
        else:
            subopt_header = 'Site_index\tOfftarget_N\tOfftarget_list\tGuide\tStar\tOligo1\tOligo2\tIsoforms\tIsoform_seqs\n'
            opt_header = 'Site_index\tGuide\tStar\tOligo1\tOligo2\tIsoforms\tIsoform_seqs\n'
        subopt.write(subopt_header)
        opt.write(opt_header)

    return subopt_name, opt_name

def add_species_to_config(config_file, species, sql_db, mrna_path):
    """
    Append a new species block to the configuration file.

    Args:
        config_file: path to config file
        specie: species name
        sql_db: SQLite database filename
        mrna_path: path to mRNA FASTA
    """

    config_path = Path(config_file) 
    # create empty file if it doesn't exist
    config_path.touch(exist_ok=True)
    block = (
        f"\n"
        f"[{species}]\n"
        f"mRNA={mrna_path.resolve()}\n"
        f"sql={sql_db.resolve()}\n"
    )

    with open(config_path, "a") as out:

        out.write(block)

    print(f"Added species '{species}' to config")

def apply_vector_to_amirna_output(data: dict, vector: str) -> None:
    """Recompute Forward/Reverse Oligo fields of a cached amiRNA output for a new vector, in place."""
    for section in ("optimal", "suboptimal", "results"):
        for entry in data.get(section, {}).values():
            fwd, rev = make_amirna_oligos(entry["amiRNA"], entry["amiRNA*"], vector)
            entry["Forward Oligo"] = fwd
            entry["Reverse Oligo"] = rev
    data["vector"] = vector


def syn_optimal_site_index(groups: dict, count: int) -> dict:
    """
    Build a {'geneset.site': guide_sequence} index of every optimal
    syn-tasiRNA site found across gene sets in a fresh pipeline run, e.g.
    {'1.1': 'ACGT...', '1.2': 'TTGC...', '2.1': 'GGCA...'}.
    """
    index = {}
    for g in range(count):
        group = groups[g]
        for o in range(1, group.get("opt", 0) + 1):
            index[f"{g + 1}.{o}"] = group["opt_r"][o]["guide"]
    return index


def syn_cached_site_index(blocks: list) -> dict:
    """
    Build the same {'geneset.site': guide_sequence} index as
    syn_optimal_site_index, but from a cached syntasirna_json 'blocks' list
    (keys there look like 'optimal 1.1' or 'result 1.1' — the label is the
    part after the last space).
    """
    index = {}
    for block in blocks:
        entries = block.get("optimal") or block.get("results") or {}
        for key, entry in entries.items():
            label = key.rsplit(" ", 1)[-1]
            index[label] = entry["syn-tasiRNA"]
    return index


def parse_syn_order(order: str, site_index: dict) -> list:
    """
    Parse a comma-separated, ordered 'geneset.site' selection (e.g.
    '1.1,3.2,2.1') into the corresponding list of guide sequences, in that
    order, validating every token against site_index.
    """
    guides = []
    for token in order.split(','):
        token = token.strip()
        if token not in site_index:
            available = ", ".join(
                sorted(site_index, key=lambda k: tuple(map(int, k.split('.'))))
            ) or "none"
            sys.exit(f"Error: site '{token}' does not exist among the optimal results found. Available: {available}.")
        guides.append(site_index[token])
    return guides


def apply_vector_to_syntasirna_output(data: dict, vector: str, target_site: str, order: str = None) -> None:
    """
    Recompute cloning_oligos of a cached syn-tasiRNA output for a new
    vector, in place, using an explicit user-chosen ordering of
    'geneset.site' selections (order, or an interactive prompt if omitted)
    instead of silently picking the first site of every gene set.
    """
    site_index = syn_cached_site_index(data.get("blocks", []))
    data["vector"] = vector
    if site_index:
        order = order or prompt_syn_order(site_index)
        syn_guides = parse_syn_order(order, site_index)
        data["selected_sites"] = [t.strip() for t in order.split(',')]
        fwd, rev = make_syntasirna_oligos(syn_guides, vector, target_site)
        data["cloning_oligos"] = {"forward": fwd, "reverse": rev}


def check_output(results_file, accession_list, output_folder, base_output=None, vector=None, target_site=None, construct=None, order=None):
    """
    If a previous full run's results already exist for this input, avoid
    re-running the (slow) TargetFinder pipeline.

    Without --vector, this keeps the original behaviour: report and exit,
    since re-running would produce the exact same output. With --vector,
    the cached canonical output (base_output) is reused and only the
    vector-specific cloning oligos are (re)computed into a new,
    vector-named output file.
    """
    if not results_file.exists():
        return

    with open(results_file, "r") as f:
        line_count = sum(1 for _ in f)

    if line_count != 4:
        return

    if not vector:
        visible_output = output_folder / f"{'_'.join(accession_list)}_psams.json"
        if base_output and base_output.exists() and not visible_output.exists():
            shutil.copy2(base_output, visible_output)
        print(f"P-SAMS already executed for {'_'.join(accession_list)}.\nResults in {output_folder.resolve()}.\nExiting script.")
        sys.exit(0)

    if not base_output.exists():
        print(f"P-SAMS already executed for {'_'.join(accession_list)}, but its cached results file ({base_output.name}) is missing. Re-running.")
        return

    with open(base_output, "r") as f:
        data = json.load(f)

    if construct == "amiRNA":
        apply_vector_to_amirna_output(data, vector)
    else:
        apply_vector_to_syntasirna_output(data, vector, target_site, order)

    vector_output = output_folder / f"{'_'.join(accession_list)}_{vector_filename_suffix(vector)}_psams.json"
    with open(vector_output, "w") as out:
        json.dump(data, out, indent=2)

    print(
        f"P-SAMS already executed for {'_'.join(accession_list)}.\n"
        f"Reusing cached results and generating cloning oligos for vector '{vector}'.\n"
        f"Output: {vector_output.resolve()}\nExiting script."
    )
    sys.exit(0)




def off_target_check(site: dict, tf_results: list, conn) -> tuple:
    """
    Use TargetFinder output to classify on-targets and off-targets.

    Parameters
    ----------
    site : dict
        Site information (must contain 'names', 'seqs', 'guide', etc.)
    tf_results : list
        Lines of output from TargetFinder
    conn : sqlite3 connection
        Database connection to annotation table

    Returns
    -------
    tuple
        (offCount, onCount, json_lines)
    """
    from html import unescape, escape
    off_target_list = []
    off_count = 0
    on_count = 0
    json_lines = []

    cursor = conn.cursor()

    for line in tf_results:

        line = line.strip()
        json_lines.append(line)

        if "Target accession" in line:
            tag, transcript = line.split(": ", 1)
            transcript = transcript.replace('"', '').replace(',', '').split()[0]

            # Query database (SQLite version)
            cursor.execute(
                "SELECT * FROM annotation WHERE transcript = ?",
                (transcript,)
            )
            result = cursor.fetchone()

            # SQLite returns tuple unless row_factory is set
            if result:
                # if using row_factory = sqlite3.Row
                try:
                    desc_value = result["description"]
                except TypeError:
                    desc_value = result[1]  # fallback: column index

                if desc_value:
                    desc = unescape(desc_value)
                    desc = escape(desc)
                    desc = desc.replace(";", "")
                    json_lines.append(f'        "Target description": "{desc}",')
                else:
                    json_lines.append('        "Target description": "unknown",')
            else:
                json_lines.append('        "Target description": "unknown",')

            # Check if on-target or off-target
            if transcript in site['names']:
                on_count += 1
            else:
                off_count += 1
                off_target_list.append(transcript)

    return off_count, on_count, json_lines, off_target_list

def base_pair(target: str, name: str, transcript: str, guide: str, construct: str):
    """
    Calculate base pairing between guide RNA and target transcript.

    Parameters
    ----------
    target : str
        Target sequence within transcript
    name : str
        Target accession or name
    transcript : str
        Full transcript sequence
    guide : str
        Guide RNA sequence
    construct : str
        Name of the construct

    Returns
    -------
    list of str
        JSON-like lines describing the hit
    """

    # Find start and end coordinates
    start = transcript.find(target)
    if start == -1:
        print(f"Warning: site {target} not found in transcript {name}!", file=sys.stderr)
        return []
    end = start + len(target) - 1

    # Base-pair scoring
    bp = {
        "AU": 0, "UA": 0, "GC": 0, "CG": 0,
        "GU": 0.5, "UG": 0.5,
        "AC": 1, "CA": 1, "AG": 1, "GA": 1,
        "UC": 1, "CU": 1,
        "A-": 1, "U-": 1, "G-": 1, "C-": 1,
        "-A": 1, "-U": 1, "-G": 1, "-C": 1,
        "AA": 1, "UU": 1, "CC": 1, "GG": 1
    }

    homology_string = ""
    cycle = 0
    score = 0
    mismatch = 0
    gu = 0

    # Convert T to U and reverse guide
    target = target.replace("T", "U")
    guide = guide.replace("T", "U")[::-1]

    guide_nts = list(guide)
    target_nts = list(target)

    while guide_nts:
        cycle += 1
        guide_base = guide_nts.pop()
        target_base = target_nts.pop()

        key = guide_base + target_base

        if cycle == 1 or cycle > 13:
            pos_score = bp.get(key, 1)
        else:
            pos_score = bp.get(key, 1) * 2

        if pos_score == 1:
            mismatch += 1
            homology_string += " "
        elif pos_score == 0.5 or pos_score == 1:
            gu += 1
            homology_string += "."
        else:
            homology_string += ":"

        score += pos_score

    homology_string = homology_string[::-1].replace(" ", "&nbsp")

    hit = [
        "      {",
        f'        "Target accession": "{name}",',
        '        "Target description": "unknown",',
        f'        "Score": "{score}",',
        f'        "Coordinates": "{start}-{end}",',
        '        "Strand": "+",',
        f'        "Target sequence": "{target}",',
        f'        "Base pairing": "{homology_string}",',
        f'        "{construct} sequence": "{guide}"',
        "      }"
    ]

    return hit

def parse_list(sep, flat_list):
    """
    Split a string by a separator or return it as a single-element list.

    Args:
        sep (str): Separator used for splitting.
        flat_list (str): Input string.

    Returns:
        list: List of parsed elements.
    """
    final_list = []

    # If separator exists in string, split it
    if sep in flat_list:
        final_list = flat_list.split(sep)
    else:
        final_list.append(flat_list)

    return final_list

def create_opt_subopt(subopt, opt, construct, no_offtarget=False, debug=False):
    """
Format and organize optimal and suboptimal target sites into structured output dictionaries, updating their labels in the TargetFinder results.

Suboptimal sites are sorted by off-target count and limited to the top three, while both groups include guide, star, oligos, and modified TF information.
Parameters
----------
subopt : list of dict
    Suboptimal sites, each containing a 'site' dictionary and an 'off_targets' count
opt : list of dict
    Optimal sites
construct : str
    Construct name used for labeling results
no_offtarget : bool, optional
    If True, off-target checking was disabled: nothing was actually
    classified as optimal, so results are labeled generically instead
debug : bool, optional
    If True, prints debugging information

Returns
-------
tuple
    opt_count : int
        Number of optimal sites processed
    subopt_count : int
        Number of suboptimal sites processed (max 3)
    opt_results : dict
        Dictionary of formatted optimal results
    subopt_results : dict
        Dictionary of formatted suboptimal results

"""
    opt_results = {}
    subopt_results = {}
    # Sort suboptimal sites by off-targets
    subopt.sort(key=lambda x: x['off_targets'])

    if debug:
        pprint(opt)

    # Process optimal results
    opt_label = "Result" if no_offtarget else "Optimal Result"
    opt_count = 0
    for i, site in enumerate(opt, start=1):
        opt_count += 1
        # Replace construct name in TF JSON
        site['tf'] = site['tf'].replace(f"{construct}{i}", f"{construct} {opt_label} {i}")

        opt_results[i] = {
            'guide': site['guide'],
            'star': site['star'],
            'oligo1': site['oligo1'],
            'oligo2': site['oligo2'],
            'tf': site['tf']
        }
  

    # Process suboptimal results
    subopt_count = 0
    for i, ssite in enumerate(subopt, start=1):
        subopt_count += 1
        site = ssite['site']
        site['tf'] = site['tf'].replace(f"{construct}{i}", f"{construct} Suboptimal Result {i}")
        subopt_results[i] = {
            'guide': site['guide'],
            'star': site['star'],
            'oligo1': site['oligo1'],
            'oligo2': site['oligo2'],
            'tf': site['tf']
        }
        if subopt_count == 3:
            break

    return opt_count, subopt_count, opt_results, subopt_results


def syntasirna_json(group_count, groups, output_file, vector=None, cloning_oligos=None, no_offtarget=False, selected_sites=None):
    """
    Build syntasiRNA JSON output from pipeline results.

    Args:
        group_count (int): number of groups
        groups (dict): pipeline result structure
        vector (str, optional): cloning vector name
        cloning_oligos (dict, optional): {"forward": ..., "reverse": ...}
        no_offtarget (bool, optional): if True, off-target checking was
            disabled, so results are labeled generically (not optimal/
            suboptimal) and the suboptimal section is omitted
        selected_sites (list, optional): ordered 'geneset.site' labels chosen
            for cloning_oligos (e.g. ['1.1', '3.2', '2.1'])
    """

    blocks = []
    set_id = 1

    for g in range(group_count):
        group = groups[g]

        section_key = "results" if no_offtarget else "optimal"
        entry_label = "result" if no_offtarget else "optimal"

        block = {"name": f"Gene set {set_id}", section_key: {}}
        if not no_offtarget:
            block["suboptimal"] = {}

        # -------------------------
        # OPTIMAL / RESULTS
        # -------------------------
        for o in range(1, group.get("opt", 0) + 1):
            entry = group["opt_r"][o]

            hits = []

            tf_raw = entry.get("tf", "")

            try:
                tf_json = json.loads(tf_raw) if isinstance(tf_raw, str) else tf_raw

                if isinstance(tf_json, dict):
                    first_key = next(iter(tf_json), None)
                    if first_key:
                        hits = tf_json[first_key].get("hits", [])

            except Exception:
                hits = []

            block[section_key][f"{entry_label} {set_id}.{o}"] = {
                "syn-tasiRNA": entry.get("guide", ""),
                "TargetFinder": hits
            }

        # -------------------------
        # SUBOPTIMAL (only when off-target checking is enabled)
        # -------------------------
        if not no_offtarget:
            for s in range(1, group.get("sub", 0) + 1):
                entry = group["sub_r"][s]

                hits = []

                tf_raw = entry.get("tf", "")

                try:
                    tf_json = json.loads(tf_raw) if isinstance(tf_raw, str) else tf_raw

                    if isinstance(tf_json, dict):
                        first_key = next(iter(tf_json), None)
                        if first_key:
                            hits = tf_json[first_key].get("hits", [])

                except Exception:
                    hits = []

                block["suboptimal"][f"suboptimal {set_id}.{s}"] = {
                    "syn-tasiRNA": entry.get("guide", ""),
                    "TargetFinder": hits
                }

        blocks.append(block)
        set_id += 1

    output = {"blocks": blocks}
    if vector:
        output["vector"] = vector
    if cloning_oligos:
        output["cloning_oligos"] = cloning_oligos
    if selected_sites:
        output["selected_sites"] = selected_sites

    with open(output_file, "w") as out:
        json.dump(output, out, indent=2)

def amirna_json(opt_count, sub_count, opt, sub, output_file, vector=None, no_offtarget=False):
    """
    Builds the JSON output for amiRNA results.

    Args:
        opt_count (int): number of optimal results
        sub_count (int): number of suboptimal results
        opt (dict): dictionary of optimal results
        sub (dict): dictionary of suboptimal results
        vector (str, optional): cloning vector name
        no_offtarget (bool, optional): if True, off-target checking was
            disabled, so results are labeled generically (not optimal/
            suboptimal) and the suboptimal section is omitted
    """
    section_key = "results" if no_offtarget else "optimal"
    entry_label = "amiRNA Result" if no_offtarget else "amiRNA Optimal Result"

    output = {section_key: {}}
    if not no_offtarget:
        output["suboptimal"] = {}
    if vector:
        output["vector"] = vector
    result_count = 0

    # Optimal / generic results
    for i in range(1, opt_count + 1):
        result_count += 1
        output[section_key][f"{entry_label} {result_count}"] = {
            "amiRNA": opt[i]["guide"],
            "amiRNA*": opt[i]["star"],
            "Forward Oligo": opt[i]["oligo1"],
            "Reverse Oligo": opt[i]["oligo2"],
            "TargetFinder": json.loads(opt[i]["tf"])
        }

    # Suboptimal results (only when off-target checking is enabled)
    if not no_offtarget:
        result_count = 0
        for i in range(1, sub_count + 1):
            try:
                tf_json = json.loads(sub[i]["tf"])
            except (json.JSONDecodeError, TypeError):
                tf_json = 'No result'  # o {}
            result_count += 1
            output["suboptimal"][f"amiRNA Suboptimal Result {result_count}"] = {
                "amiRNA": sub[i]["guide"],
                "amiRNA*": sub[i]["star"],
                "Forward Oligo": sub[i]["oligo1"],
                "Reverse Oligo": sub[i]["oligo2"],
               "TargetFinder": tf_json
            }

    with open(output_file, "w") as out:
        json.dump(output, out, indent=2)

