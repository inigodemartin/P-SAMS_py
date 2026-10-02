#!/usr/bin/env python3
"""
Simple syn-tasiRNA oligo designer

This script generates the sense and antisense oligos to clone one or more
syn-tasiRNAs following the format:

Sense oligo:     TTTA + syn-tasiRNA_1 + syn-tasiRNA_2 + ...
Antisense oligo: CCGA + reverse_complement(last) + ... + reverse_complement(first)

Based on the cloning logic described in Supplemental Protocol 1 for AtTAS1c-D2-B/c vectors.
"""

from typing import List, Tuple


DNA_COMPLEMENT = str.maketrans("ATGCatgc", "TACGtacg")


def clean_sequence(seq: str) -> str:
    """Normalize sequence: remove spaces and convert U to T."""
    seq = seq.strip().replace(" ", "").replace("-", "").upper().replace("U", "T")
    return seq


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def validate_syn_tasirna(seq: str) -> None:
    """Validate that the syn-tasiRNA sequence is a 21-nt DNA sequence."""
    if len(seq) != 21:
        raise ValueError(f"Sequence '{seq}' must be 21 nt long. Current length: {len(seq)}")
    invalid = set(seq) - set("ATGC")
    if invalid:
        raise ValueError(f"Sequence '{seq}' contains invalid bases: {', '.join(sorted(invalid))}")


def build_oligos(entries: List[Tuple[str, str]]) -> Tuple[str, str]:
    """
    Build oligos from a list of (name, sequence) tuples.

    Sense oligo:
        TTTA + seq1 + seq2 + ...

    Antisense oligo:
        CCGA + revcomp(seqN) + revcomp(seqN-1) + ...
    """
    sequences = [seq for _, seq in entries]
    sense = "TTTA" + "".join(sequences)
    antisense = "CCGA" + "".join(reverse_complement(seq) for seq in reversed(sequences))
    return sense, antisense


def get_user_entries() -> List[Tuple[str, str]]:
    """Interactively collect syn-tasiRNA names and sequences."""
    print("syn-tasiRNA oligo designer")
    print("-" * 30)
    print("Enter one or more syn-tasiRNAs.")
    print("Type an empty name to finish.\n")

    entries: List[Tuple[str, str]] = []

    while True:
        name = input("syn-tasiRNA name: ").strip()
        if not name:
            break

        seq = clean_sequence(input(f"Sequence for {name} (21 nt): "))
        validate_syn_tasirna(seq)
        entries.append((name, seq))
        print()

    if not entries:
        raise ValueError("No syn-tasiRNAs were entered.")

    return entries


def main() -> None:
    try:
        entries = get_user_entries()
        sense, antisense = build_oligos(entries)

        print("\nRESULTS")
        print("=" * 30)
        print("Inserted syn-tasiRNAs:")
        for i, (name, seq) in enumerate(entries, start=1):
            print(f"{i}. {name}: {seq}")

        print("\nSense oligo:")
        print(sense)
        print(f"Length: {len(sense)} nt")

        print("\nAntisense oligo:")
        print(antisense)
        print(f"Length: {len(antisense)} nt")

    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()
