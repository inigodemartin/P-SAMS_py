from src.utils import design_guide_RNA, oligo_designer, base_pair, _status_loop, off_target_check
import sqlite3
import json
import time
import threading
import subprocess
from concurrent.futures import ThreadPoolExecutor

def get_tsites(ids: dict, seed: int, bg=False, conn=None, debug=False) -> list:
    """
    Identify all putative target sites in transcripts.

    Parameters
    ----------
    ids : dict
        {transcript_id: sequence}
    seed : int
        Length of the seed region
    bg : bool
        If True, filter k-mers present in background database
    conn : mysql connection, optional
        Required if bg=True
    debug : bool
        Print debug info

    Returns
    -------
    list of dict
        Each dict contains {'name': transcript_id, 'seq': 21-nt site}
    """
    not_t_sites = []
    t_sites = []
    site_length = 21
    offset = site_length - seed - 1
    if debug:
        print("Finding sites in foreground transcripts...")

    discard = set()  # kmers to discard
    found = set()    # kmers already checked

    cursor = None
    if bg and conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

    for transcript, seq in ids.items():
        length = len(seq)
        if debug:
            print(f"Transcript {transcript} is {length} nt long")

        for i in range(length - site_length + 1):
            site = seq[i:i + site_length]
            if len(site) < site_length:
                continue
            kmer = site[offset:offset + seed]
            # Background filtering
            ## En la base de datos tenemos todos los posibles kmers (de 15 nt) del genoma a analizar y los transcritos donde aparecen
            ## En esta parte del codigo se chequea que los kmers de nuestras secuencias únicamente estén en nuestro input
            if bg and cursor:
                is_bg = False

                if kmer in discard:
                    is_bg = True
                    continue

                elif kmer not in found:
                    cursor.execute("SELECT * FROM kmers WHERE kmer = ?", (kmer,))
                    results = cursor.fetchall()

                    for result in results:
                        accessions = result["transcripts"].split(',')# aqui sacamos los transcritos que contienen ese k-mer

                        for accession in accessions:
                            if accession not in ids: # Si ese id no se encuentra en nuestros ids de input se descarta (osea, sería un offtarget)
                                is_bg = True
                                discard.add(kmer) # Aqui se almacenan los kmers que tienen off-targets
                                break
                    found.add(kmer) #Aqui se almacenan los kmers que ya se han chequeado
                    if is_bg:
                        continue         
            # Store site
            t_sites.append({'name': transcript, 'seq': site})
    return t_sites


def group_tsites(t_sites: list, seed: int, debug: bool = False) -> list:
    """
    Group target sites based on seed sequence.

    Parameters
    ----------
    t_sites : list of dict
        [{'name': transcript_id, 'seq': sequence}, ...]
    seed : int
        Length of the seed region
    debug : bool
        Print debug information

    Returns
    -------
    list of dict
        [{'names': 'id1;id2', 'seqs': 'seq1;seq2'}, ...]
    """

    site_length = 21
    offset = site_length - seed - 1

    if debug:
        print("Grouping sites...")

    # --------------------------------------------------
    # Sort by seed sequence, then by transcript name
    # --------------------------------------------------
    t_sites = sorted(
        t_sites,
        key=lambda x: (x['seq'][offset:offset + seed], x['name'])
    )

    gsites = []
    names = []
    seqs = []
    last_seq = None
    
    for i, row in enumerate(t_sites):
        current_seq = row['seq']
        current_seed = current_seq[offset:offset + seed]

        # Initialize first group
        if not names:
            last_seq = current_seq

        last_seed = last_seq[offset:offset + seed]

        if current_seed == last_seed:
            names.append(row['name'])
            seqs.append(current_seq)

        else:
            # Save previous group
            gsites.append({
                'names': ";".join(names),
                'seqs': ";".join(seqs)
            })

            # Start new group
            names = [row['name']]
            seqs = [current_seq]
            last_seq = current_seq

    # --------------------------------------------------
    # Add last group 
    # --------------------------------------------------
    if names:
        gsites.append({
            'names': ";".join(names),
            'seqs': ";".join(seqs)
        })

    if debug:
        print("done")
    

    # LO QUE HACE ESTA FUNCION ES AGRUPAR TODAS LAS TARGETSITES QUE TIENEN LA MISMA SEED
    return gsites



def score_sites(gsites: list, min_site_count: int, seed: int, fb: str, debug: bool = False) -> list:
    """
    Score grouped target sites based on sequence similarity and positional rules.

    Parameters
    ----------
    gsites : list of dict
        [{'names': ..., 'seqs': ...}, ...]
    min_site_count : int
        Minimum number of sequences per group
    seed : int
        Seed length
    fb : str
        Foldback type
    debug : bool
        Debug flag

    Returns
    -------
    list of dict
        Scored sites
    """

    site_length = 21
    offset = site_length - seed - 1
    scored = []

    for site in gsites:

        sites = site['seqs'].split(';')

        if len(sites) < min_site_count:

            continue
        # -----------------------------
        # Position 21 (index 20)
        # -----------------------------
        nts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        
        for seq in sites:
            nts[seq[20]] += 1
        if nts['A'] == len(sites):
            site['p21'] = 1
        elif nts['G'] == len(sites):
            site['p21'] = 2
        elif nts['A'] + nts['G'] == len(sites):
            site['p21'] = 3
        else:
            site['p21'] = 4

        # -----------------------------
        # Position 1 (index 0)
        # -----------------------------
        nts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        for seq in sites:
            nts[seq[0]] += 1
        if any(nts[n] == len(sites) for n in nts):
            site['p1'] = 1
        elif nts['G'] > 0 and nts['T'] > 0:
            site['p1'] = 2
        else:
            site['p1'] = 1
        # -----------------------------
        # Position 2 (index 1)
        # -----------------------------
        nts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        for seq in sites:
            nts[seq[1]] += 1

        if any(nts[n] == len(sites) for n in nts):
            site['p2'] = 1
        elif nts['G'] > 0 and nts['A'] > 0 and nts['G'] + nts['A'] == len(sites):
            site['p2'] = 2
        else:
            site['p2'] = 3
        # -----------------------------
        # Position 3 (index 2)
        # -----------------------------
        nts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        for seq in sites:
            nts[seq[2]] += 1

        if nts['G'] == len(sites):
            site['p3'] = 1
        elif any(nts[n] == len(sites) for n in ['A', 'C', 'T']):
            site['p3'] = 2
        else:
            site['p3'] = 3

        # -----------------------------
        # Other mismatches
        # -----------------------------
        site['other_mm'] = 0

        for i in range(3, offset):
            nts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
            for seq in sites:
                nts[seq[i]] += 1

            if not any(nts[n] == len(sites) for n in nts):
                site['other_mm'] += 1

        # -----------------------------
        # Guide + oligos
        # -----------------------------
        guide = design_guide_RNA(site)
       
        if len(guide) < site_length:
            continue

        site['guide'] = guide
        
        star, oligo1, oligo2 = oligo_designer(guide, fb)

        site['star'] = star
        site['oligo1'] = oligo1
        site['oligo2'] = oligo2

        scored.append(site)

    if debug:
        print("Sorting and outputing results...")

    # -----------------------------
    # Sorting
    # -----------------------------
    scored = sorted(
        scored,
        key=lambda x: (
            x['other_mm'],
            x['p21'],
            x['p3'],
            x['p2'],
            x['p1']
        )
    )

    if debug:
        print(f"Analyzing {len(scored)} total sites...")

    return scored



def _run_targetfinder(targetfinder, guide, mRNA_fa, query_name):
    """
    Run TargetFinder for a single guide and return its stdout as a list of lines.
    """
    cmd = [
        targetfinder,
        "-s", guide,
        "-d", mRNA_fa,
        "-q", query_name,
        "-p", "json"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.splitlines()


def serial_jobs(target_count, construct, ids, site_scores,targetfinder, mRNA_fa, conn, bg,subopt_name, opt_name, output_folder,accession_list, tf_dir, potential_target_n,unlimit=False, fasta=False, jobs=1):
    """
    Run serial jobs for amiRNA evaluation.

    Parameters
    ----------
    target_count : int
        Expected number of on-targets.

    construct : str
        Construct prefix name.

    ids : dict
        Dictionary containing transcript sequences:
        {transcript_id: sequence}

    bg : bool
        Whether to check off-targets using TargetFinder.

    site_scores : list
        List of grouped and scored target sites.

    targetfinder : str
        Path to TargetFinder script.

    mRNA_fa : str
        Path to mRNA database.

    fasta : bool
        Whether input comes from FASTA.

    unlimit : bool
        If True, do not limit the number of optimal results.

    jobs : int
        Number of TargetFinder calls to run in parallel. Default = 1 (serial).

    Returns
    -------
    tuple
        (opt, subopt)
    """
    with open(subopt_name, 'a', buffering=1) as subopt_out, open(opt_name, 'a', buffering=1) as opt_out:

        start_time = time.time()

        optimal_ref = [0]
        suboptimal_ref = [0]

        stop_event = threading.Event()

        thread = threading.Thread(
            target=_status_loop,
            args=(start_time, accession_list, optimal_ref, suboptimal_ref, stop_event, potential_target_n),
            daemon=True
        )
        thread.start()

        opt = []
        subopt = []
        count = 0

        n_sites = len(site_scores)
        jobs = max(1, jobs)

        # ------------------------------------------------------------
        # When jobs > 1, run TargetFinder for several candidates ahead
        # of time in a thread pool (subprocess calls release the GIL),
        # while still consuming results in order. Each prefetched call
        # uses a placeholder query name (since the final "{construct}N"
        # name depends on how many optimal sites have been found so
        # far, which is only known once results are processed in
        # order); the placeholder is swapped for the real name once
        # the result is consumed.
        # ------------------------------------------------------------
        executor = None
        pending = {}

        if bg and jobs > 1:
            executor = ThreadPoolExecutor(max_workers=jobs)

            def _submit(i):
                placeholder = f"{construct}_job{i}"
                pending[i] = (
                    placeholder,
                    executor.submit(_run_targetfinder, targetfinder, site_scores[i]['guide'], mRNA_fa, placeholder)
                )

            for i in range(min(jobs, n_sites)):
                _submit(i)

        next_submit = jobs

        for idx, site in enumerate(site_scores):

            count += 1
            site['name'] = f"{construct}{optimal_ref[0]}"

            if bg:

                if executor:
                    placeholder, future = pending.pop(idx)
                    tf_results = future.result()

                    # Keep the pipeline full: queue the next candidate as
                    # soon as a slot frees up.
                    if next_submit < n_sites:
                        _submit(next_submit)
                        next_submit += 1

                    tf_results = [line.replace(placeholder, site['name']) for line in tf_results]
                else:
                    tf_results = _run_targetfinder(targetfinder, site['guide'], mRNA_fa, site['name'])

                # TargetFinder prints "No results for <query>" (not JSON) when the
                # designed guide doesn't align even to its own target above the
                # score cutoff. Such candidates are unusable, skip them.
                if not tf_results or tf_results[0].strip().startswith("No results for"):
                    continue

                off_targets, on_targets, json_data, offtarget_list = off_target_check(
                    site, tf_results, conn
                )

                site['offtarget_list'] = offtarget_list
                site['offtarget_number'] = len(offtarget_list)

                if fasta:

                    insert = []
                    seqs = site['seqs'].split(';')
                    names = site['names'].split(';')

                    for seq, name in zip(seqs, names):
                        hit = base_pair(seq, name, ids[name], site['guide'], construct)
                        insert.append("\n      ".join(hit))

                    if off_targets == 0:

                        site['tf'] = "\n".join([
                            '{',
                            f'  "{construct}{optimal_ref[0]}": {{',
                            '    "hits": [',
                            ",\n".join(insert),
                            '    ]',
                            '  }',
                            '}'
                        ])

                        opt.append(site)

                        opt_out.write(
                            f"{count}\t{site['guide']}\t{site['star']}\t"
                            f"{site['oligo1']}\t{site['oligo2']}\t"
                            f"{site['names']}\t{site['seqs']}\n"
                        )

                        with open(f"{tf_dir}/site_{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        optimal_ref[0] += 1

                    else:

                        new_json = json_data[:3]
                        new_json.append(",\n".join(insert) + ',')
                        new_json.extend(json_data[3:])

                        site['tf'] = new_json

                        subopt.append({
                            'off_targets': off_targets,
                            'site': site
                        })

                        subopt_out.write(
                            f"{count}\t{site['offtarget_number']}\t{site['offtarget_list']}\t"
                            f"{site['guide']}\t{site['star']}\t"
                            f"{site['oligo1']}\t{site['oligo2']}\t"
                            f"{site['names']}\t{site['seqs']}\n"
                        )

                        with open(f"{tf_dir}/site_{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        suboptimal_ref[0] += 1

                else:

                    site['tf'] = "\n".join(json_data)

                    if off_targets == 0 and on_targets == target_count:

                        opt.append(site)

                        opt_out.write(
                            f"{count}\t{site['guide']}\t{site['star']}\t"
                            f"{site['oligo1']}\t{site['oligo2']}\t"
                            f"{site['names']}\t{site['seqs']}\n"
                        )

                        with open(f"{tf_dir}/site_{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        optimal_ref[0] += 1

                    else:

                        subopt.append({
                            'off_targets': off_targets,
                            'site': site
                        })

                        subopt_out.write(
                            f"{count}\t{site['offtarget_number']}\t{site['offtarget_list']}\t"
                            f"{site['guide']}\t{site['star']}\t"
                            f"{site['oligo1']}\t{site['oligo2']}\t"
                            f"{site['names']}\t{site['seqs']}\n"
                        )

                        with open(f"{tf_dir}/site_{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        suboptimal_ref[0] += 1

            else:

                insert = []
                seqs = site['seqs'].split(';')
                names = site['names'].split(';')

                for seq, name in zip(seqs, names):
                    hit = base_pair(seq, name, ids[name], site['guide'], construct)
                    insert.append("\n      ".join(hit))

                site['tf'] = "\n".join([
                    '{',
                    f'  "{site["name"]}": {{',
                    '    "hits": [',
                    ",\n      ".join(insert),
                    '    ]',
                    '  }',
                    '}'
                ])

                opt.append(site)
                optimal_ref[0] += 1

            if optimal_ref[0] == 3 and not unlimit:
                break

        if executor:
            executor.shutdown(wait=False, cancel_futures=True)

        stop_event.set()
        thread.join()

    print()
    return opt, subopt