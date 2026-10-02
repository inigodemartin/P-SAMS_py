#!/usr/bin/env python3
"""
Minimal syn-tasiRNA oligo designer for pENTR-B/c or pMDC32-B/c.

The user provides:
1. A 22-nt miRNA target site sequence.
2. One or more 21-nt syn-tasiRNA sequences.

Oligo architecture:

Forward:
TGTA + 22-nt target site + TAGACCATTTA + syn-tasiRNA_1 + syn-tasiRNA_2 + ...

Reverse:
AATG + reverse_complement(syn-tasiRNA_N) + ... + reverse_complement(syn-tasiRNA_1)
     + TAAATGGTCTA + reverse_complement(22-nt target site)

Note:
For two syn-tasiRNAs, the final oligos are 79 nt, as in the example.
For one or more than two syn-tasiRNAs, the length changes accordingly.
"""

from typing import List, Tuple


DNA_COMPLEMENT = str.maketrans("ATGCatgc", "TACGtacg")
FORWARD_PREFIX = "TGTA"
REVERSE_PREFIX = "AATG"
SPACER_FORWARD = "TAGACCATTTA"
SPACER_REVERSE = "TAAATGGTCTA"


def clean_sequence(seq: str) -> str:
    """Normalize sequence: remove spaces/hyphens, convert U to T, uppercase."""
    return seq.strip().replace(" ", "").replace("-", "").upper().replace("U", "T")


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def validate_dna(seq: str, expected_length: int, label: str) -> None:
    """Validate sequence length and DNA alphabet."""
    if len(seq) != expected_length:
        raise ValueError(f"{label} must be {expected_length} nt long. Current length: {len(seq)} nt.")
    invalid = set(seq) - set("ATGC")
    if invalid:
        raise ValueError(f"{label} contains invalid bases: {', '.join(sorted(invalid))}")


def build_oligos(target_site: str, syn_tasirnas: List[Tuple[str, str]]) -> Tuple[str, str]:
    """
    Build forward and reverse oligos.

    Forward:
    TGTA + target_site + TAGACCATTTA + syn-tasiRNAs

    Reverse:
    AATG + reverse-complemented syn-tasiRNAs in reverse order
         + TAAATGGTCTA + reverse-complement target_site
    """
    syn_sequences = [seq for _, seq in syn_tasirnas]

    forward = (
        FORWARD_PREFIX
        + target_site
        + SPACER_FORWARD
        + "".join(syn_sequences)
    )

    reverse = (
        REVERSE_PREFIX
        + "".join(reverse_complement(seq) for seq in reversed(syn_sequences))
        + SPACER_REVERSE
        + reverse_complement(target_site)
    )

    return forward, reverse


def get_user_input() -> Tuple[str, str, List[Tuple[str, str]]]:
    """Collect vector, target site, and syn-tasiRNA entries from the user."""
    print("Minimal syn-tasiRNA oligo designer")
    print("-" * 38)
    print("Compatible vectors: pENTR-B/c or pMDC32-B/c\n")

    vector = input("Vector name [pMDC32-B/c]: ").strip() or "pMDC32-B/c"
    allowed_vectors = {"pENTR-B/c", "pMDC32-B/c"}
    if vector not in allowed_vectors:
        raise ValueError(f"Vector must be one of: {', '.join(sorted(allowed_vectors))}")

    target_site = clean_sequence(input("22-nt miRNA target site sequence: "))
    validate_dna(target_site, 22, "Target site")

    print("\nEnter one or more syn-tasiRNAs.")
    print("Type an empty name to finish.\n")

    syn_tasirnas: List[Tuple[str, str]] = []

    while True:
        name = input("syn-tasiRNA name: ").strip()
        if not name:
            break

        seq = clean_sequence(input(f"Sequence for {name} (21 nt): "))
        validate_dna(seq, 21, name)
        syn_tasirnas.append((name, seq))
        print()

    if not syn_tasirnas:
        raise ValueError("At least one syn-tasiRNA is required.")

    return vector, target_site, syn_tasirnas


def main() -> None:
    try:
        vector, target_site, syn_tasirnas = get_user_input()
        forward, reverse = build_oligos(target_site, syn_tasirnas)

        construct_name = "_".join(name for name, _ in syn_tasirnas)
        forward_name = f"{construct_name}_F"
        reverse_name = f"{construct_name}_R"

        print("\nRESULTS")
        print("=" * 38)
        print(f"Vector: {vector}")
        print(f"Target site: {target_site}")
        print(f"Spacer sequence: {SPACER_FORWARD}")
        print(f"Number of syn-tasiRNAs: {len(syn_tasirnas)}\n")

        for i, (name, seq) in enumerate(syn_tasirnas, start=1):
            print(f"{i}. {name}: {seq}")

        print("\nForward oligonucleotide:")
        print(f"Name: {forward_name}")
        print(f"Sequence: {forward}")
        print(f"Length: {len(forward)} nt")

        print("\nReverse oligonucleotide:")
        print(f"Name: {reverse_name}")
        print(f"Sequence: {reverse}")
        print(f"Length: {len(reverse)} nt")

    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()
