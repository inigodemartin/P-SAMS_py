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
