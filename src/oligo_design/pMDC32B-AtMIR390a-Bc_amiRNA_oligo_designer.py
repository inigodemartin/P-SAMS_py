#!/usr/bin/env python3

def normalize_seq(seq: str) -> str:
    return seq.upper().replace(" ", "").replace("U", "T")

def revcomp(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]

def validate_dna(seq: str, name: str) -> None:
    allowed = set("ACGT")
    bad = set(seq) - allowed
    if bad:
        raise ValueError(f"{name} contains invalid bases: {''.join(sorted(bad))}")

def make_oligos(amirna: str, amirna_star: str):
    amirna = normalize_seq(amirna)
    amirna_star = normalize_seq(amirna_star)

    if len(amirna) != 21:
        raise ValueError(f"amiRNA must be 21 nt, got {len(amirna)} nt")
    if len(amirna_star) != 21:
        raise ValueError(f"amiRNA* must be 21 nt, got {len(amirna_star)} nt")

    validate_dna(amirna, "amiRNA")
    validate_dna(amirna_star, "amiRNA*")

    upstream_star = revcomp(amirna[-2:])
    star_block = upstream_star + amirna_star[:19]

    forward = "TGTA" + amirna + "ATGATGATCACATTCGTTATCTATTTTTT" + star_block
    reverse = "AATG" + revcomp(forward[4:])

    return forward, reverse, upstream_star, star_block

if __name__ == "__main__":
    name = input("Enter amiRNA name: ").strip()
    amirna = input("Enter amiRNA (21 nt): ").strip()
    amirna_star = input("Enter amiRNA* (21 nt): ").strip()

    forward, reverse, upstream_star, star_block = make_oligos(amirna, amirna_star)

    amirna = normalize_seq(amirna)
    amirna_star = normalize_seq(amirna_star)

    amirna_star_name = f"{name}*"
    forward_name = f"{name}_F"
    reverse_name = f"{name}_R"

    print(f"\namiRNA name: {name}")
    print(f"amiRNA sequence: {amirna}")
    print(f"amiRNA* name: {amirna_star_name}")
    print(f"amiRNA* sequence: {amirna_star}")
    print(f"2 variable nt upstream of amiRNA*: {upstream_star}")
    print(f"Star block used in oligos: {star_block}")
    print(f"Forward oligo name: {forward_name}")
    print(f"Forward oligo sequence: {forward}")
    print(f"Reverse oligo name: {reverse_name}")
    print(f"Reverse oligo sequence: {reverse}")

    filename = f"{name}_oligos_output.txt"
    with open(filename, "w") as f:
        f.write(f"amiRNA name: {name}\n")
        f.write(f"amiRNA sequence: {amirna}\n\n")
        f.write(f"amiRNA* name: {amirna_star_name}\n")
        f.write(f"amiRNA* sequence: {amirna_star}\n\n")
        f.write(f"2 variable nt upstream of amiRNA*: {upstream_star}\n")
        f.write(f"Star block used in oligos: {star_block}\n\n")
        f.write(f"Forward oligo name: {forward_name}\n")
        f.write(f"Forward oligo sequence: {forward}\n\n")
        f.write(f"Reverse oligo name: {reverse_name}\n")
        f.write(f"Reverse oligo sequence: {reverse}\n")

    print(f"\nResults saved to: {filename}")
