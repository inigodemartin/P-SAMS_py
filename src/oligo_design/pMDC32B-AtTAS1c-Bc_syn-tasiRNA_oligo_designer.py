#!/usr/bin/env python3

"""
syn-tasiRNA oligo generator for pMDC32B-AtTAS1c-B/c

Rules:
- Forward: ATTA + concatenated syn-tasiRNA sequences
- Reverse: GTTC + reverse complement of the concatenated sequence
- Forward and reverse must have identical length
"""

FWD_OVERHANG = "ATTA"
REV_OVERHANG = "GTTC"
VECTOR = "pMDC32B-AtTAS1c-B/c"

def clean(seq):
    return seq.upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("-", "")

def revcomp(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]

def validate_syn(seq, name="syn-tasiRNA"):
    if len(seq) != 21:
        raise ValueError(f"{name} must be 21 nt, got {len(seq)} nt.")
    bad = sorted(set(seq) - set("ACGT"))
    if bad:
        raise ValueError(f"{name} contains invalid base(s): {''.join(bad)}")

def generate_oligos(syn_tasirnas):
    syn_tasirnas = [clean(s) for s in syn_tasirnas]
    if not syn_tasirnas:
        raise ValueError("Provide at least one syn-tasiRNA.")

    for i, s in enumerate(syn_tasirnas, 1):
        validate_syn(s, f"syn-tasiRNA-{i}")

    core = "".join(syn_tasirnas)

    forward = FWD_OVERHANG + core
    reverse = REV_OVERHANG + revcomp(core)

    if len(forward) != len(reverse):
        raise RuntimeError("Forward and reverse oligos must have the same length.")

    return {
        "vector": VECTOR,
        "core_sequence": core,
        "forward_oligo": forward,
        "reverse_oligo": reverse,
        "length_nt": len(forward),
    }

if __name__ == "__main__":
    # Example
    syns = [
        "TCTTGTAACGCGCTTTCCCAG",
        "TTCGCTGTACAGTTCTTTCGC",
    ]
    res = generate_oligos(syns)
    print("Vector:", res["vector"])
    print("Forward:", res["forward_oligo"])
    print("Reverse:", res["reverse_oligo"])
    print("Length:", res["length_nt"], "nt")
