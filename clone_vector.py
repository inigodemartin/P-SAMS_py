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


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate vector-specific cloning oligos for an existing psams.py run."
    )
    ap.add_argument("-o", "--output_folder", required=True,
        help="The '..._psams_output' folder from a previous psams.py run.")
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

    args = ap.parse_args()

    if args.target_site:
        ts = args.target_site.strip().upper().replace("U", "T")
        invalid = set(ts) - set("ACGT")
        if invalid or len(ts) != 22:
            ap.error("--target-site must be a 22-nt DNA sequence (only A, C, G, T bases allowed).")
        args.target_site = ts

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
    output_folder = Path(args.output_folder).resolve()
    if not output_folder.exists():
        sys.exit(f"Error: output folder not found: {output_folder}")

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
