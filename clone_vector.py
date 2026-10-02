#!/usr/bin/env python3

"""
Generate vector-specific cloning oligos for a P-SAMS run that already found
its optimal sites (see psams.py). Doesn't rerun TargetFinder — it only
(re)computes the cloning oligos for the given vector and, for syntasiRNA,
the chosen site order.

This is a separate step on purpose: for syntasiRNA, which sites to clone
and in what order is a decision only the user can make, and can only be
made once the optimal sites are already known — so it can't happen before
or during the (possibly long) pipeline run in psams.py.
"""

import argparse
import json
import sys
from pathlib import Path

from src.utils import syn_cached_site_index, apply_vector_to_amirna_output, apply_vector_to_syntasirna_output, write_amirna_tsv, write_syntasirna_tsv
from src.oligo_design import AMIRNA_VECTORS, SYNTASIRNA_VECTORS, prompt_target_site, vector_filename_suffix
from src.insert_design import (INSERT_VECTORS, amirna_cached_site_index, build_insert,
                               check_foldback, design_basename, parse_modules,
                               spec_uses_cache, write_exports)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate vector-specific cloning oligos for an existing psams.py run."
    )
    ap.add_argument("-o", "--output_folder",
        help="The '..._psams_output' folder from a previous psams.py run. For an insert "
             "vector whose modules are all typed out, this is just where the files are "
             "written, and defaults to the current directory.")
    ap.add_argument("-V", "--vector", required=True,
        help="Cloning vector name (see psams.py's README for the list per construct).")
    ap.add_argument("-O", "--order",
        metavar="SITES",
        help="syntasiRNA only. Ordered, comma-separated list of sites to clone into the vector, "
             "as 'geneset.site' (e.g. '1.1,3.2,2.1'). If omitted, you're prompted interactively.")
    ap.add_argument("-T", "--target-site",
        dest="target_site",
        metavar="SEQ",
        help="syntasiRNA only. 22-nt miRNA target site sequence, required for the pMDC32B-B/c "
             "vector. If omitted and needed, it will be requested interactively.")
    ap.add_argument("-M", "--modules",
        metavar="SPEC",
        help="Insert vectors only (" + ", ".join(INSERT_VECTORS) + "). Ordered, comma-separated "
             "list of the modules to lay down. Each sequence is either addressed into this "
             "run's optimal results or typed out in full, so the vectors can be designed with "
             "no run behind them. Tokens: 'ami:N' (cached amiRNA result N); "
             "'ami:<name>:<21-nt amiRNA>[:<21-nt amiRNA*>]' (typed out; the amiRNA* is derived "
             "from the amiRNA when left out); 'mir173a[:<name>]' (the fixed AtMIR173a "
             "precursor); 'syn:<target site>:<guides>' (a syn-tasiRNA module, target site named "
             "e.g. NbmiR482a, spelled out as 22 nt, or 'name=<22 nt>'; guides '+'-joined, each "
             "a cached 'geneset.site' reference, a 'name=<21 nt>' pair, or bare 21 nt). "
             "Examples: 'ami:1,syn:NbmiR482a:1.1+2.1' or "
             "'ami:amiR-NbSu:TGTATGACTCCCGGAATTCCA'.")
    ap.add_argument("-n", "--name",
        help="Insert vectors only. Construct name, used in the exports and the output filenames. "
             "Defaults to the run's accession key, or to the vector's own default when there is "
             "no run behind the design.")

    args = ap.parse_args()

    if args.target_site:
        ts = args.target_site.strip().upper().replace("U", "T")
        invalid = set(ts) - set("ACGT")
        if invalid or len(ts) != 22:
            ap.error("--target-site must be a 22-nt DNA sequence (only A, C, G, T bases allowed).")
        args.target_site = ts

    is_insert = args.vector in INSERT_VECTORS
    if is_insert and not args.modules:
        ap.error(f"--modules is required for the insert vector '{args.vector}'.")
    if not is_insert and not args.output_folder:
        ap.error("--output_folder is required.")
    if args.modules and not is_insert:
        ap.error(
            f"--modules only applies to the insert vectors ({', '.join(INSERT_VECTORS)}); "
            f"'{args.vector}' is cloned from an oligo pair instead."
        )

    return args


def _find_cache_json(output_folder: Path, vector: str) -> Path:
    cache_dir = output_folder / ".cache"
    matches = sorted(cache_dir.glob("*_psams.json")) if cache_dir.exists() else []
    if not matches:
        sys.exit(f"Error: no cached results found in {cache_dir} — run psams.py on this input first.")

    if len(matches) == 1:
        return matches[0]

    # amiRNA and syntasiRNA can now be run into the same output folder
    # (psams.py namespaces their cache files separately), so more than one
    # cached run may be sitting here. Vector names are unique to one
    # construct, so the requested -V vector disambiguates which one to use.
    wants_syntasirna = vector in SYNTASIRNA_VECTORS
    wants_amirna = vector in AMIRNA_VECTORS
    for path in matches:
        with open(path) as f:
            is_syn = "blocks" in json.load(f)
        if is_syn and wants_syntasirna:
            return path
        if not is_syn and wants_amirna:
            return path

    sys.exit(
        f"Error: multiple cached runs found in {cache_dir} "
        f"({', '.join(p.name for p in matches)}) and none match vector '{vector}'."
    )


def _load_insert_caches(output_folder: Path) -> tuple:
    """
    Load both of a run's cached results for an insert vector.

    An insert can combine amiRNA and syn-tasiRNA modules in one construct,
    and those come from two separate psams.py runs with a cache file each,
    so unlike _find_cache_json there is nothing to disambiguate here: read
    whichever of the two are present and let the module specification
    decide what it needs.
    """
    cache_dir = output_folder / ".cache"
    matches = sorted(cache_dir.glob("*_psams.json")) if cache_dir.exists() else []
    if not matches:
        sys.exit(f"Error: no cached results found in {cache_dir} — run psams.py on this input first.")

    amirna_data = None
    syn_data = None
    accession_keys = []

    for path in matches:
        with open(path) as f:
            data = json.load(f)

        run_key = path.name[: -len("_psams.json")]
        for tag in ("_syntasiRNA", "_amiRNA"):
            if run_key.endswith(tag):
                run_key = run_key[: -len(tag)]
                break
        accession_keys.append(run_key)

        if "blocks" in data:
            syn_data = data
        else:
            amirna_data = data

    # Only the eudicot AtMIR390a foldback has a defined insert-vector design,
    # so a monocot run must not be wrapped in these flanks (see
    # insert_design.check_foldback). Runs cached before the field existed
    # carry no "foldback" key and are taken at the default, eudicot.
    check_foldback((amirna_data or {}).get("foldback"))

    return amirna_data, syn_data, accession_keys[0]


def _syntasirna_vector_output_path(output_folder: Path, accession_key: str, vector: str, new_oligos: dict) -> Path:
    """
    Pick where to write this run's syntasiRNA cloning oligos for `vector`.

    The same gene set(s) can be cloned more than once with a different site
    selection/order (-O), which produces different oligos each time — so a
    second, different result must never silently overwrite the first one.
    The first ever output for a given vector keeps the plain, unsuffixed
    name; subsequent *different* oligo sets get a numeric suffix (_1, _2,
    ...). Re-running with the exact same selection/order (same oligos)
    reuses/overwrites its own matching file instead of piling up duplicates.
    """
    suffix = vector_filename_suffix(vector)
    n = 0
    while True:
        name = (
            f"{accession_key}_{suffix}_psams.json" if n == 0
            else f"{accession_key}_{suffix}_psams_{n}.json"
        )
        path = output_folder / name
        if not path.exists():
            return path
        with open(path) as f:
            existing = json.load(f)
        if existing.get("cloning_oligos") == new_oligos:
            return path
        n += 1


def main():
    args = parse_args()
    output_folder = Path(args.output_folder or ".").resolve()
    if not output_folder.exists():
        sys.exit(f"Error: output folder not found: {output_folder}")

    if args.vector in INSERT_VECTORS:
        # A specification that addresses nothing into a previous run needs no
        # run at all: the tool then stands on its own and only needs somewhere
        # to write.
        if spec_uses_cache(args.modules):
            amirna_data, syn_data, accession_key = _load_insert_caches(output_folder)
            amirna_index = amirna_cached_site_index(amirna_data or {})
            syn_index = syn_cached_site_index((syn_data or {}).get("blocks", []))
        else:
            amirna_index, syn_index, accession_key = {}, {}, None

        modules = parse_modules(args.modules, amirna_index, syn_index)
        design = build_insert(args.vector, modules, args.name or accession_key or "")

        base = design_basename(design)
        json_path = output_folder / f"{base}_design.json"
        with open(json_path, "w") as out:
            json.dump(design, out, indent=2)
        written = write_exports(design, output_folder)

        print(
            f"Design generated successfully ({len(modules)} module"
            f"{'s' if len(modules) > 1 else ''}, insert {design['insert_length']} bp)."
        )
        for path in [json_path] + written:
            print(f"Output: {path}")
        return

    cache_json = _find_cache_json(output_folder, args.vector)
    with open(cache_json) as f:
        data = json.load(f)

    accession_key = cache_json.name[: -len("_psams.json")]
    is_syntasirna = "blocks" in data

    vectors = SYNTASIRNA_VECTORS if is_syntasirna else AMIRNA_VECTORS
    if args.vector not in vectors:
        available = ", ".join(vectors)
        sys.exit(f"Error: '{args.vector}' is not a valid vector for this construct. Available: {available}.")

    if is_syntasirna:
        site_index = syn_cached_site_index(data.get("blocks", []))
        if not site_index:
            sys.exit("Error: no optimal syn-tasiRNA sites were found for this run; nothing to clone.")

        target_site = args.target_site
        if args.vector == "pMDC32B-B/c" and not target_site:
            target_site = prompt_target_site()

        apply_vector_to_syntasirna_output(data, args.vector, target_site, args.order)
        vector_output = _syntasirna_vector_output_path(output_folder, accession_key, args.vector, data["cloning_oligos"])
    else:
        apply_vector_to_amirna_output(data, args.vector)
        vector_output = output_folder / f"{accession_key}_{vector_filename_suffix(args.vector)}_psams.json"

    with open(vector_output, "w") as out:
        json.dump(data, out, indent=2)

    if is_syntasirna:
        write_syntasirna_tsv(data, vector_output.with_suffix(".tsv"))
    else:
        write_amirna_tsv(data, vector_output.with_suffix(".tsv"))

    print(f"Cloning oligos generated for vector '{args.vector}'.\nOutput: {vector_output}")


if __name__ == "__main__":
    main()
