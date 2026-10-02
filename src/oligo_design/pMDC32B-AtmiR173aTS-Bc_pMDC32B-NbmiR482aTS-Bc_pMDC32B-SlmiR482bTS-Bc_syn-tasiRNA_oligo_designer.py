#!/usr/bin/env python3

VECTORS = {
    "AtmiR173a": "pMDC32-AtmiR173aTS-B/c",
    "NbmiR482a": "pMDC32-NbmiR482aTS-B/c",
    "SlmiR482b": "pMDC32-SlmiR482bTS-B/c",
}

FWD_OVERHANG = "TTTA"
REV_OVERHANG = "CCGA"

def clean(seq):
    """Normalize DNA/RNA sequence."""
    return seq.upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("-", "")

def revcomp(seq):
    """Reverse complement of a DNA sequence."""
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]

def validate_syn_tasirna(seq, name="syn-tasiRNA"):
    if len(seq) != 21:
        raise ValueError(f"{name} must be 21 nt, got {len(seq)} nt.")
    invalid = sorted(set(seq) - set("ACGT"))
    if invalid:
        raise ValueError(f"{name} contains invalid base(s): {''.join(invalid)}")

def generate_oligos(syn_tasirnas, target_site):
    """
    Generate oligos for one or more syn-tasiRNAs.

    Forward oligo:
        TTTA + concatenated syn-tasiRNA sequence(s)

    Reverse oligo:
        CCGA + reverse complement of the concatenated syn-tasiRNA sequence(s)

    Therefore, forward and reverse oligos always have the same length.
    """

    if target_site not in VECTORS:
        raise ValueError(f"Target site must be one of: {', '.join(VECTORS)}")

    syn_tasirnas = [clean(s) for s in syn_tasirnas]

    if not syn_tasirnas:
        raise ValueError("At least one syn-tasiRNA sequence is required.")

    for i, seq in enumerate(syn_tasirnas, 1):
        validate_syn_tasirna(seq, f"syn-tasiRNA-{i}")

    core = "".join(syn_tasirnas)

    forward = FWD_OVERHANG + core
    reverse = REV_OVERHANG + revcomp(core)

    if len(forward) != len(reverse):
        raise RuntimeError("Forward and reverse oligos do not have the same length.")

    return {
        "vector": VECTORS[target_site],
        "core_sequence": core,
        "forward_oligo": forward,
        "reverse_oligo": reverse,
        "forward_length": len(forward),
        "reverse_length": len(reverse),
    }


if __name__ == "__main__":
    # Example
    syn_tasirnas = [
        "TCTTGTAACGCGCTTTCCCAG",
        "TTCGCTGTACAGTTCTTTCGC",
    ]

    result = generate_oligos(syn_tasirnas, "NbmiR482a")

    print("Selected vector:", result["vector"])
    print("Forward oligo:", result["forward_oligo"])
    print("Reverse oligo:", result["reverse_oligo"])
    print("Forward length:", result["forward_length"], "nt")
    print("Reverse length:", result["reverse_length"], "nt")
