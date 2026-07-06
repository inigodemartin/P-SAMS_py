def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _make_oligos(fwd_prefix: str, rev_prefix: str, body: str):
    return fwd_prefix + body, rev_prefix + _revcomp(body)


# (fwd_prefix, rev_prefix, spacer)
AMIRNA_VECTORS = {
    "pMDC32B-AtMIR390a-B/c":        ("TGTA", "AATG", "ATGATGATCACATTCGTTATCTATTTTTT"),
    "pMDC32B-OsMIR390-B/c":          ("CTTG", "CATG", "ATGATGATCACATTCGTTATCTATTTTTT"),
    "pMDC32B-BS-AtMIR390a-B/c":      ("TGTA", "AATG", "CGAAATCAAACT"),
    "pMDC32B-BS-AtMIR390a-A18G-B/c": ("TGTG", "AATG", "CGAAATCAAACT"),
}

# (fwd_prefix, rev_prefix, target_site_spacer or None)
# When target_site_spacer is not None, a 22-nt miRNA target site must be provided.
SYNTASIRNA_VECTORS = {
    "pMDC32B-B/c":               ("TGTA", "AATG", "TAGACCATTTA"),
    "pMDC32B-AtTAS1c-B/c":       ("ATTA", "GTTC", None),
    "pMDC32B-AtTAS1c-D2-B/c":    ("TTTA", "CCGA", None),
    "pMDC32B-AtmiR173aTS-B/c":   ("TTTA", "CCGA", None),
    "pMDC32B-NbmiR482aTS-B/c":   ("TTTA", "CCGA", None),
    "pMDC32B-SlmiR482bTS-B/c":   ("TTTA", "CCGA", None),
}


def make_amirna_oligos(amirna: str, amirna_star: str, vector: str):
    """Return (forward_oligo, reverse_oligo) for an amiRNA candidate and cloning vector."""
    fwd_prefix, rev_prefix, spacer = AMIRNA_VECTORS[vector]
    upstream_star = _revcomp(amirna[-2:])
    star_block = upstream_star + amirna_star[:19]
    body = amirna + spacer + star_block
    return _make_oligos(fwd_prefix, rev_prefix, body)


def make_syntasirna_oligos(syn_guides: list, vector: str, target_site: str = None):
    """Return (forward_oligo, reverse_oligo) for a list of syn-tasiRNA sequences and cloning vector."""
    fwd_prefix, rev_prefix, ts_spacer = SYNTASIRNA_VECTORS[vector]
    if ts_spacer is not None:
        body = target_site + ts_spacer + "".join(syn_guides)
    else:
        body = "".join(syn_guides)
    return _make_oligos(fwd_prefix, rev_prefix, body)


def select_vector(construct: str) -> str:
    """Interactively prompt the user to pick a cloning vector for the given construct type."""
    vectors = list(AMIRNA_VECTORS if construct == "amiRNA" else SYNTASIRNA_VECTORS)

    print("\nAvailable cloning vectors:")
    for i, name in enumerate(vectors, start=1):
        print(f"  {i}. {name}")

    while True:
        choice = input(f"Select a vector [1-{len(vectors)}]: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(vectors):
            return vectors[int(choice) - 1]
        print(f"Invalid selection '{choice}'. Please enter a number between 1 and {len(vectors)}.")


def prompt_target_site() -> str:
    """Interactively prompt the user for the 22-nt miRNA target site sequence."""
    while True:
        seq = input("Enter the 22-nt miRNA target site sequence: ").strip().upper().replace("U", "T")
        if len(seq) == 22 and set(seq) <= set("ACGT"):
            return seq
        print("Invalid sequence: must be a 22-nt DNA sequence (A, C, G, T only).")


def vector_filename_suffix(vector: str) -> str:
    """Return a filesystem-safe suffix for a vector name (e.g. 'pMDC32B-B/c' -> 'pMDC32B-Bc')."""
    return vector.replace("/", "")
