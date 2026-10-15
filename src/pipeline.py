from src.utils import design_guide_RNA, oligo_designer, base_pair, _status_loop, _render_status_line, off_target_check
import sqlite3
import json
import time
import threading
import subprocess
import contextlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

def get_tsites(ids: dict, seed: int, bg=False, backgrounds=None, debug=False) -> list:
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
    backgrounds : list of utils.Background, optional
        Every background transcriptome to screen against. Required if
        bg=True. A seed found in a non-target transcript of ANY of them is
        discarded, so this first filter is already the intersection of all
        the hosts a guide must be clean in.
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

    cursors = []
    if bg and backgrounds:
        for background in backgrounds:
            background.conn.row_factory = sqlite3.Row
            cursors.append(background.conn.cursor())

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
            if bg and cursors:
                is_bg = False

                if kmer in discard:
                    is_bg = True
                    continue

                elif kmer not in found:
                    # Every background is queried, not just the first: the
                    # seed has to be absent from all of them to survive.
                    for cursor in cursors:
                        cursor.execute("SELECT * FROM kmers WHERE kmer = ?", (kmer,))
                        results = cursor.fetchall()

                        for result in results:
                            accessions = result["transcripts"].split(',')# aqui sacamos los transcritos que contienen ese k-mer

                            for accession in accessions:
                                if accession not in ids: # Si ese id no se encuentra en nuestros ids de input se descarta (osea, sería un offtarget)
                                    is_bg = True
                                    discard.add(kmer) # Aqui se almacenan los kmers que tienen off-targets
                                    break
                            if is_bg:
                                break
                        if is_bg:
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



def _load_tf_cache(cache_path: Path) -> dict:
    """
    Load the shared, construct-agnostic TargetFinder result cache from a
    previous run against the same gene(s)/species in this output folder
    (see serial_jobs). TargetFinder's result for a given guide only
    depends on that guide and the background transcriptome, never on
    which construct designed it — so an amiRNA run and a syntasiRNA run
    for the same gene(s) can reuse each other's TargetFinder calls instead
    of re-invoking the (slow) external tool for guides already tested.
    """
    results = {}
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                raw = json.load(f)
            # Entries without a "background" come from a cache written
            # before runs could screen against several transcriptomes.
            # There is no record of which transcriptome produced them, so
            # they are dropped rather than risk reusing a result computed
            # against a different background.
            results = {
                _tf_cache_key(entry["background"], guide_key.rsplit("\t", 1)[-1]): (entry["query_name"], entry["lines"])
                for guide_key, entry in raw.items()
                if isinstance(entry, dict) and "background" in entry
            }
        except Exception:
            results = {}
    return {"path": cache_path, "lock": threading.Lock(), "results": results}


def _tf_cache_key(background_name: str, guide: str) -> str:
    """
    A TargetFinder result depends on the guide AND on the transcriptome it
    was searched against, so both make up the cache key. Keying on the
    guide alone (as earlier versions did, when the background was always
    the same single species) would silently hand back a result computed
    against another host.
    """
    return f"{background_name}\t{guide}"


def _tf_cache_store(tf_cache: dict, background_name: str, guide: str, query_name: str, lines: list) -> None:
    with tf_cache["lock"]:
        tf_cache["results"][_tf_cache_key(background_name, guide)] = (query_name, lines)
        try:
            raw = {}
            for key, (qn, ln) in tf_cache["results"].items():
                bg_name, guide_seq = key.split("\t", 1)
                raw[key] = {"background": bg_name, "guide": guide_seq, "query_name": qn, "lines": ln}
            with open(tf_cache["path"], "w") as f:
                json.dump(raw, f)
        except Exception:
            pass  # best-effort cache; a write failure shouldn't break the run


def _run_targetfinder(targetfinder, guide, background, query_name, tf_cache=None):
    """
    Run TargetFinder for a single guide against one background
    transcriptome and return its stdout as a list of lines. If tf_cache is
    given (see _load_tf_cache) and already has a result for this exact
    guide AND background, that cached result is relabeled to query_name and
    returned instead of re-invoking TargetFinder — the same label-swap
    already used for prefetched threaded results below.
    """
    if tf_cache is not None:
        cached = tf_cache["results"].get(_tf_cache_key(background.name, guide))
        if cached is not None:
            cached_query_name, lines = cached
            return [line.replace(cached_query_name, query_name) for line in lines]

    cmd = [
        targetfinder,
        "-s", guide,
        "-d", background.mRNA_fa,
        "-q", query_name,
        "-p", "json"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.splitlines()

    if tf_cache is not None:
        _tf_cache_store(tf_cache, background.name, guide, query_name, lines)

    return lines


def _no_results(lines: list) -> bool:
    """
    TargetFinder prints "No results for <query>" (not JSON) when nothing in
    the searched transcriptome aligns to the guide above the score cutoff.
    """
    return not lines or lines[0].strip().startswith("No results for")


def _hits_from_lines(lines: list) -> list:
    """
    Turn one TargetFinder JSON document (as a list of lines, possibly with
    annotation descriptions already injected by off_target_check) into its
    flat list of hit objects.

    Results from several backgrounds have to be merged into a single
    document, and each invocation emits a complete document of its own, so
    they are parsed and recombined as objects rather than spliced as text.
    """
    try:
        document = json.loads("\n".join(lines))
    except ValueError:
        return []
    if not isinstance(document, dict):
        return []
    first_key = next(iter(document), None)
    if first_key is None:
        return []
    return document[first_key].get("hits", [])


def _tag_species(hit: dict, background_name: str) -> dict:
    """
    Record which transcriptome a hit came from, right next to the accession
    it refers to. Only used when screening against more than one, so a
    single-species run's output is unchanged.
    """
    tagged = {}
    for key, value in hit.items():
        tagged[key] = value
        if key == "Target accession":
            tagged["Target species"] = background_name
    return tagged


def _tf_document(name: str, hits: list) -> str:
    """Wrap hits back into the {query_name: {"hits": [...]}} document the rest of the pipeline stores."""
    return json.dumps({name: {"hits": hits}}, indent=4)


def serial_jobs(target_count, construct, ids, site_scores,targetfinder, backgrounds, bg,subopt_name, opt_name, output_folder,accession_list, tf_dir, potential_target_n,unlimit=False, fasta=False, jobs=1, limit=3, gene_set=None, existing_opt=None, existing_subopt=None, seen_guides=None, start_count=0):
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

    backgrounds : list of utils.Background
        Every background transcriptome to screen candidates against. A
        candidate is only optimal when it has no off-target in any of them,
        which is what lets one guide be designed against a pathogen gene
        and checked against all of its hosts at once.

    fasta : bool
        Whether input comes from FASTA.

    unlimit : bool
        If True, do not limit the number of optimal results.

    limit : int
        Maximum number of optimal results to find before stopping.
        Ignored when unlimit is True. Default = 3.

    gene_set : int, optional
        1-indexed gene set number for syntasiRNA multi-gene-set runs. When
        set, every TSV row is prefixed with a Gene_set column so rows from
        different gene sets (which share the same TSV file) can be told
        apart. None (default) for amiRNA, which has no gene sets.

    existing_opt, existing_subopt : list, optional
        Already-evaluated optimal/suboptimal candidates from a previous
        partial run (see utils.load_resume_state), to resume from instead
        of starting over.

    seen_guides : set, optional
        Guide sequences already evaluated in a previous partial run;
        candidates with these guides are skipped instead of being sent to
        TargetFinder again.

    start_count : int, optional
        Highest Site_index/tf_results file number already used by a
        previous partial run, so resumed rows/files continue numbering
        instead of colliding with existing ones.

    jobs : int
        Number of TargetFinder calls to run in parallel. Default = 1 (serial).

    Returns
    -------
    tuple
        (opt, subopt)
    """
    # In no-offtarget mode (bg=False), every candidate is treated as a plain
    # "result" (see the bg-False branch below) — nothing is ever classified
    # as optimal/suboptimal, so the TSV files are never written to. Skip
    # creating them at all rather than leaving behind header-only files.
    subopt_cm = open(subopt_name, 'a', buffering=1) if bg else contextlib.nullcontext(None)
    opt_cm = open(opt_name, 'a', buffering=1) if bg else contextlib.nullcontext(None)
    row_prefix = f"{gene_set}\t" if gene_set is not None else ""
    gs_prefix = f"gs{gene_set}_" if gene_set is not None else ""

    def _opt_row(site, count):
        # syntasiRNA (gene_set is not None) checkpoint rows never carry
        # Oligo1/Oligo2 — see create_outputs.
        cols = [f"{row_prefix}{count}", site['guide'], site['star']]
        if gene_set is None:
            cols += [site['oligo1'], site['oligo2']]
        cols += [site['names'], site['seqs']]
        return "\t".join(cols) + "\n"

    def _subopt_row(site, count):
        cols = [
            f"{row_prefix}{count}", str(site['offtarget_number']), str(site['offtarget_list']),
            site['guide'], site['star'],
        ]
        if gene_set is None:
            cols += [site['oligo1'], site['oligo2']]
        cols += [site['names'], site['seqs']]
        return "\t".join(cols) + "\n"

    # Resuming: candidates whose guide was already evaluated in a previous
    # partial run are dropped from site_scores so TargetFinder never runs
    # on them again; only genuinely new candidates go through the loop
    # below, and their results are appended to the already-found ones.
    if seen_guides:
        site_scores = [s for s in site_scores if s['guide'] not in seen_guides]

    with subopt_cm as subopt_out, opt_cm as opt_out:

        start_time = time.time()

        opt = list(existing_opt) if existing_opt else []
        subopt = list(existing_subopt) if existing_subopt else []

        # Already have enough optimal results from a previous run: nothing
        # left to do for this gene set.
        if not unlimit and len(opt) >= limit:
            site_scores = []

        optimal_ref = [len(opt)]
        suboptimal_ref = [len(subopt)]

        stop_event = threading.Event()

        thread = threading.Thread(
            target=_status_loop,
            args=(start_time, accession_list, optimal_ref, suboptimal_ref, stop_event, potential_target_n, bg),
            daemon=True
        )
        thread.start()

        count = start_count

        n_sites = len(site_scores)
        jobs = max(1, jobs)

        # Shared across constructs (see _load_tf_cache): a guide already
        # TargetFinder-tested by e.g. an earlier amiRNA run for this same
        # gene(s) doesn't need to be tested again for a later syntasiRNA
        # run in the same output folder, and vice versa.
        tf_cache = _load_tf_cache(Path(output_folder) / ".cache" / "tf_results_cache.json") if bg else None

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

        backgrounds = backgrounds or []
        # Screening against several transcriptomes is the only case that
        # needs each hit labelled with the host it was found in; a
        # single-background run keeps exactly the output it had before.
        multi = len(backgrounds) > 1
        offtarget_hosts = Counter()

        lookahead = jobs

        if bg and jobs > 1:
            executor = ThreadPoolExecutor(max_workers=jobs)

            def _submit(i):
                placeholder = f"{construct}_job{i}"
                pending[i] = (
                    placeholder,
                    [
                        executor.submit(_run_targetfinder, targetfinder, site_scores[i]['guide'],
                                        background, placeholder, tf_cache)
                        for background in backgrounds
                    ]
                )

            # Every candidate now costs one TargetFinder call per
            # background, so fewer candidates are kept in flight to still
            # run -j/--jobs external calls at a time.
            lookahead = max(1, jobs // max(1, len(backgrounds)))
            for i in range(min(lookahead, n_sites)):
                _submit(i)

        next_submit = lookahead

        for idx, site in enumerate(site_scores):

            count += 1
            site['name'] = f"{construct}{optimal_ref[0]}"

            if bg:

                if executor:
                    placeholder, futures = pending.pop(idx)
                    per_background = [
                        (background, [line.replace(placeholder, site['name']) for line in future.result()])
                        for background, future in zip(backgrounds, futures)
                    ]

                    # Keep the pipeline full: queue the next candidate as
                    # soon as a slot frees up.
                    if next_submit < n_sites:
                        _submit(next_submit)
                        next_submit += 1

                else:
                    per_background = [
                        (background, _run_targetfinder(targetfinder, site['guide'], background, site['name'], tf_cache))
                        for background in backgrounds
                    ]

                searched = [(background, lines) for background, lines in per_background if not _no_results(lines)]

                # A background returning nothing means different things in
                # the two input modes. With -a the target is itself one of
                # the searched transcripts, so a guide that aligns nowhere
                # does not even hit its own target and is unusable. With -f
                # the target may well be absent from every transcriptome
                # searched (a pathogen gene screened against its hosts), so
                # finding nothing there is the best possible outcome, and
                # the on-target pairing is computed separately below.
                if not searched and not fasta:
                    continue

                off_targets = 0
                on_targets = 0
                offtarget_list = []
                tf_hits = []

                for background, lines in searched:
                    off, on, json_data, off_list = off_target_check(
                        site, lines, background.conn, construct
                    )
                    hits = _hits_from_lines(json_data)

                    if multi:
                        off_list = [f"{background.name}:{accession}" for accession in off_list]
                        hits = [_tag_species(hit, background.name) for hit in hits]

                    off_targets += off
                    on_targets += on
                    offtarget_list.extend(off_list)
                    tf_hits.extend(hits)

                    if off:
                        offtarget_hosts[background.name] += 1

                site['offtarget_list'] = offtarget_list
                site['offtarget_number'] = len(offtarget_list)

                if fasta:

                    on_hits = []
                    seqs = site['seqs'].split(';')
                    names = site['names'].split(';')

                    for seq, name in zip(seqs, names):
                        hit = base_pair(seq, name, ids[name], site['guide'], construct)
                        # base_pair returns nothing when the site cannot be
                        # located in its own transcript; such a hit has no
                        # JSON to merge.
                        if hit:
                            on_hits.append(json.loads("\n".join(hit)))

                    # Same as before: the searched transcriptomes' own hits
                    # are only carried over when there is something
                    # off-target to show, since the on-target pairing is
                    # already described by on_hits.
                    site['tf'] = _tf_document(site['name'], on_hits + (tf_hits if off_targets else []))

                    if off_targets == 0:

                        opt.append(site)

                        opt_out.write(_opt_row(site, count))

                        with open(f"{tf_dir}/site_{gs_prefix}{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        optimal_ref[0] += 1

                    else:

                        subopt.append({
                            'off_targets': off_targets,
                            'site': site
                        })

                        subopt_out.write(_subopt_row(site, count))

                        with open(f"{tf_dir}/site_{gs_prefix}{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        suboptimal_ref[0] += 1

                else:

                    site['tf'] = _tf_document(site['name'], tf_hits)

                    if off_targets == 0 and on_targets == target_count:

                        opt.append(site)

                        opt_out.write(_opt_row(site, count))

                        with open(f"{tf_dir}/site_{gs_prefix}{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        optimal_ref[0] += 1

                    else:

                        subopt.append({
                            'off_targets': off_targets,
                            'site': site
                        })

                        subopt_out.write(_subopt_row(site, count))

                        with open(f"{tf_dir}/site_{gs_prefix}{count:04d}_TargetFinder_result.json", "w") as siteout:
                            json.dump(json.loads(site['tf']), siteout, indent=4)

                        suboptimal_ref[0] += 1

            else:

                insert = []
                seqs = site['seqs'].split(';')
                names = site['names'].split(';')

                for seq, name in zip(seqs, names):
                    hit = base_pair(seq, name, ids[name], site['guide'], construct)
                    if hit:
                        insert.append(json.loads("\n".join(hit)))

                site['tf'] = _tf_document(site['name'], insert)

                opt.append(site)
                optimal_ref[0] += 1

            if optimal_ref[0] == limit and not unlimit:
                break

        if executor:
            # wait=True: any candidates already prefetched and running when
            # the limit was hit must finish before this function returns —
            # with wait=False those TargetFinder subprocesses kept running
            # as orphaned, non-daemon executor threads, which Python's
            # interpreter shutdown then silently blocked on (the process
            # looked hung, or dumped a KeyboardInterrupt traceback, well
            # after "Finished running P-SAMS successfully!" had printed).
            executor.shutdown(wait=True, cancel_futures=True)

        stop_event.set()
        thread.join()

        # The for loop above can drain a whole burst of already-completed
        # futures (jobs > 1) faster than the once-a-second status thread
        # can observe, so the last line it printed may be stale (e.g. still
        # showing 0 optimal sites right before this gene set's true count).
        # Render one final, synchronous line with the actual end-of-loop
        # counts before moving on.
        _render_status_line(start_time, accession_list, optimal_ref, suboptimal_ref, potential_target_n, bg)

    print()

    # With several hosts, knowing WHICH one rejected a candidate is what
    # tells the user whether dropping a host would recover usable guides.
    if multi and offtarget_hosts:
        print("  Candidates with off-targets per background:")
        for background in backgrounds:
            print(f"    {background.name}: {offtarget_hosts[background.name]}")

    return opt, subopt