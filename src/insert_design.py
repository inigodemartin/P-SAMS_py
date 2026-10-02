"""
Build full inserts for the vectors that are ordered as a synthetic fragment
instead of being cloned from a pair of annealed oligos.

The vectors in src/oligo_design.py all share one shape: a single art-sRNA
goes into a pMDC32B-*-B/c backbone, and the output is a short oligo pair to
anneal. The vectors here are the other shape: one or more modules are laid
end to end between two fixed flanks, and the output is the whole insert,
meant to be ordered as a synthetic fragment.

Three vectors use it:

  PVX, TRV                  viral vectors, flanked by plain 20-nt homology
                            arms with no restriction site
  pMDC32B-B/c-multimodule   the B/c backbone entered through BsaI Golden
                            Gate, which is what makes more than one module
                            per construct possible

A module is one of:

  amiRNA        an AtMIR390a-style foldback carrying an amiRNA/amiRNA* duplex
  AtMIR173a     a fixed 117-nt foldback, inserted verbatim
  syn-tasiRNA   a 22-nt miRNA target site followed by one or more 21-nt
                syn-tasiRNAs

Modules are butt-joined: nothing is inserted between them. The only internal
linker is the AtTAS1c-derived spacer sitting inside a syn-tasiRNA module,
between its target site and its syn-tasiRNAs.

Sequences come either from a previous psams.py run's cached optimal results
or, for standalone use with no run behind them, straight from the command
line. Everything the Oligo Designer Suite HTML does for these three vectors
is reproduced here: the same fixed sequences, the same amiRNA* derivation,
the same validation, and byte-identical TXT/FASTA/CSV exports.
"""

import csv
import io
import re
import sys


# ---------------------------------------------------------------------------
# Fixed sequences
# ---------------------------------------------------------------------------

# (5' flank, 3' flank)
INSERT_VECTORS = {
    # Viral vectors: plain homology/overlap arms, no restriction site.
    "PVX": ("AGAGGTCAGCACCAGCTAGC", "AGGGTTTGTTAAGTTTCCCT"),
    "TRV": ("CACTTACCCGAGTTAACGCC", "ATGTCCCGAAGACATTAAAC"),
    # GCATCG + GGTCTC + A + TGTA  /  CATT + A + GAGACC + CGATGC
    # GGTCTC and its reverse complement GAGACC are the BsaI sites; TGTA and
    # CATT are the overhangs BsaI releases.
    "pMDC32B-B/c-multimodule": ("GCATCGGGTCTCATGTA", "CATTAGAGACCCGATGC"),
}

# AtMIR390a foldback, around the amiRNA/amiRNA* duplex.
AMI_BS5 = "AGTAGAGAAGAATCTGTA"   # 5' basal stem
AMI_DSL = "CGAAATCAAACT"         # distal stem-loop
AMI_BS3 = "TTGGCTCTTCTTACT"      # 3' basal stem

AMI_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
# Transition (A<->G, C<->T), used for the single central bulge in the star.
AMI_TRANS = {"G": "A", "C": "T", "A": "G", "T": "C"}

# AtTAS1c-derived spacer, between a target site and its syn-tasiRNAs.
SYN_SPACER = "TAGACCATTTA"

# Fixed AtMIR173a precursor, used verbatim as a whole module.
AT_MIR173A = (
    "TGGTGATTAAGTACTTTCGCTTGCAGAGAGAAATCACAGTGGTCAAAAAAGTTGTAGTTTT"
    "CTTAAAGTCTCTTTCCTCTGTGATTCTCTGTGTAAGCGAAAGAGCTTGCTCCCTAA"
)

# 22-nt miRNA target sites that can be named instead of typed out. The value
# is (label as the HTML writes it, sequence).
BUILT_IN_TARGET_SITES = {
    "AtmiR173a": ("AtmiR173a TS", "GTGATTTTTCTCTACAAGCGAA"),
    "NbmiR482a": ("NbmiR482a TS", "GTGGTATGGGGGGAGTCGGGAA"),
    "SlmiR482b": ("SlmiR482b TS", "GGCATGGGCGGTGTAGGCAAGA"),
}

# Only the eudicot AtMIR390a foldback has a defined insert-vector design.
# psams.py can also build monocot (OsMIR390) foldbacks, whose basal stems and
# star differ; wrapping one in the constants above would silently produce a
# chimeric precursor, so it is refused rather than guessed at.
SUPPORTED_FOLDBACK = "eudicot"


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------

def clean_sequence(seq: str, label: str) -> str:
    """
    Normalise a hand-typed sequence: drop whitespace and dashes, uppercase,
    read RNA as DNA, and refuse anything that is not A/C/G/T/U.
    """
    cleaned = re.sub(r"[\s-]", "", seq.strip()).upper().replace("U", "T")
    if not re.fullmatch(r"[ACGT]*", cleaned):
        sys.exit(f"Error: invalid sequence for {label}: {seq}. Use only A/C/G/T/U.")
    return cleaned


def _require_length(seq: str, length: int, label: str) -> str:
    if len(seq) != length:
        sys.exit(f"Error: {label} must be {length} nt. Received {len(seq)}.")
    return seq


def ami_star(amirna: str) -> str:
    """
    Derive the 21-nt amiRNA* from the amiRNA: 19 nt complementary to it, with
    the central bulge at star position 9 set by transition from amiRNA
    position 11, plus the fixed AtMIR390a 3' overhang "CA".
    """
    _require_length(amirna, 21, "amiRNA")
    bases = [
        AMI_TRANS[amirna[10]] if k == 9 else AMI_COMP[amirna[19 - k]]
        for k in range(1, 20)
    ]
    return "".join(bases) + "CA"


def check_foldback(foldback) -> None:
    """Refuse a cached run whose amiRNAs are not eudicot AtMIR390a ones."""
    if foldback and foldback != SUPPORTED_FOLDBACK:
        sys.exit(
            f"Error: this run was designed with the '{foldback}' foldback, and the "
            f"insert vectors are only defined for the '{SUPPORTED_FOLDBACK}' "
            f"AtMIR390a foldback. Its basal stems and amiRNA* differ, so the "
            f"insert would not fold as intended."
        )


# ---------------------------------------------------------------------------
# Modules
# ---------------------------------------------------------------------------

def amirna_module(name: str, amirna: str, amirna_star: str = None) -> dict:
    """
    One AtMIR390a-style foldback carrying an amiRNA/amiRNA* duplex.

    amirna_star is the 21-nt star as the rest of P-SAMS reports it: 19 nt
    facing the amiRNA plus the fixed "CA" 3' overhang. Left out, it is
    derived from the amiRNA. The two nucleotides that sit between the distal
    stem-loop and the star are scaffold, not part of the star, so they are
    derived here from the amiRNA's own 3' dinucleotide rather than being
    carried around.
    """
    amirna = _require_length(amirna, 21, "amiRNA")
    if amirna_star is None:
        amirna_star = ami_star(amirna)
    else:
        _require_length(amirna_star, 21, "amiRNA*")

    seq = AMI_BS5 + amirna + AMI_DSL + _revcomp(amirna[-2:]) + amirna_star + AMI_BS3
    return {
        "type": "amiRNA",
        "name": name,
        "amiRNA": amirna,
        "amiRNA*": amirna_star,
        "sequence": seq,
    }


def mir173a_module(name: str = "AtMIR173a") -> dict:
    """The fixed AtMIR173a precursor, inserted verbatim."""
    return {"type": "AtMIR173a", "name": name, "sequence": AT_MIR173A}


def syntasirna_module(target_name: str, target_site: str, guides: list) -> dict:
    """
    One syn-tasiRNA module: a 22-nt miRNA target site, the AtTAS1c-derived
    spacer, then every 21-nt syn-tasiRNA end to end.

    guides is a list of (name, sequence) pairs, in the order they are to be
    laid down.
    """
    target_site = _require_length(target_site, 22, "target site")
    guides = list(guides)
    if not guides:
        sys.exit("Error: a syn-tasiRNA module requires at least one syn-tasiRNA.")
    for i, (name, seq) in enumerate(guides, start=1):
        _require_length(seq, 21, f"syn-tasiRNA {i}")

    seq = target_site + SYN_SPACER + "".join(s for _, s in guides)
    return {
        "type": "syn-tasiRNA",
        "target_name": target_name,
        "target_site": target_site,
        "guides": guides,
        "sequence": seq,
    }


# ---------------------------------------------------------------------------
# Insert
# ---------------------------------------------------------------------------

# shape -> (Type: label, CSV type token, filename stem fallback, filename tag)
_DESIGN_LABELS = {
    ("PVX", "ami"): (
        "PVX amiRNA (shc / AtMIR390a-based)",
        "PVX_amiRNA", "pvx_amiRNA", "PVX_amiRNA"),
    ("PVX", "syn"): (
        "PVX syn-tasiRNA (AtTAS1c-based)",
        "PVX_syn_tasiRNA", "pvx_syn_tasiRNA", "PVX_syn_tasiRNA"),
    ("PVX", "hybrid"): (
        "PVX hybrid (amiRNA + syn-tasiRNA in tandem)",
        "PVX_hybrid_amiRNA_syn_tasiRNA", "pvx_hybrid", "PVX_hybrid"),
    ("TRV", "ami"): (
        "TRV amiRNA (shc / MIR390-based)",
        "TRV_amiRNA", "trv_amiRNA", "TRV_amiRNA"),
    ("TRV", "syn"): (
        "TRV syn-tasiRNA (AtTAS1c-based)",
        "TRV_syn_tasiRNA", "trv_syn_tasiRNA", "TRV_syn_tasiRNA"),
    ("TRV", "hybrid"): (
        "TRV Multimodule (amiRNA + syn-tasiRNA in tandem)",
        "TRV_multimodule_amiRNA_syn_tasiRNA", "trv_multimodule", "TRV_multimodule"),
}

# The B/c multimodule backbone is one panel in the HTML whatever it carries,
# so its labels do not vary with the module mix.
_BC_LABELS = (
    "B/c multimodule (amiRNA + syn-tasiRNA in tandem)",
    "Bc_multimodule_amiRNA_syn_tasiRNA", "bc_multimodule", "Bc_multimodule")


def design_shape(modules: list) -> str:
    """
    Which of the HTML's three insert panels this module mix corresponds to.

    AtMIR173a counts as a multimodule construct rather than an amiRNA one:
    it is a MIR173 precursor, not a MIR390 foldback, and the HTML only ever
    offers it inside the hybrid/multimodule panels.
    """
    kinds = {m["type"] for m in modules}
    has_syn = "syn-tasiRNA" in kinds
    has_ami = "amiRNA" in kinds
    has_mir173 = "AtMIR173a" in kinds

    if has_syn and (has_ami or has_mir173):
        return "hybrid"
    if has_syn:
        return "syn"
    if has_mir173:
        return "hybrid"
    return "ami"


def _labels(vector: str, modules: list) -> tuple:
    if vector.startswith("pMDC32B"):
        return _BC_LABELS
    return _DESIGN_LABELS[(vector, design_shape(modules))]


def build_insert(vector: str, modules: list, construct_name: str) -> dict:
    """
    Lay every module between the vector's flanks and return the finished
    design: the insert itself plus everything the exports need to describe
    it.
    """
    if vector not in INSERT_VECTORS:
        raise KeyError(vector)
    if not modules:
        sys.exit("Error: an insert needs at least one module.")

    flank5, flank3 = INSERT_VECTORS[vector]
    insert = (flank5 + "".join(m["sequence"] for m in modules) + flank3).upper()
    type_label, csv_type, fallback, tag = _labels(vector, modules)

    return {
        "construct_name": construct_name,
        "vector": vector,
        "type": type_label,
        "csv_type": csv_type,
        "name_fallback": fallback,
        "file_tag": tag,
        # The B/c insert is ordered as a single-stranded fragment; the viral
        # ones as double-stranded. Only the wording of the report differs.
        "strandedness": "ssDNA" if vector.startswith("pMDC32B") else "dsDNA",
        "flank_5": flank5.upper(),
        "flank_3": flank3.upper(),
        "modules": modules,
        "insert": insert,
        "insert_length": len(insert),
    }


# ---------------------------------------------------------------------------
# Module specification
# ---------------------------------------------------------------------------

def amirna_cached_site_index(data: dict) -> dict:
    """
    Build a {index: {'amiRNA': ..., 'amiRNA*': ...}} map of the optimal
    amiRNA results in a cached amiRNA run, keyed by the result number its
    label ends with (see utils.amirna_json, which writes labels such as
    'amiRNA Optimal Result 1').
    """
    index = {}
    for section in ("optimal", "results"):
        for label, entry in (data.get(section) or {}).items():
            match = re.search(r"(\d+)\s*$", label)
            if match:
                index[int(match.group(1))] = entry
    return index


def _resolve_target_site(token: str) -> tuple:
    """
    A target site is named from BUILT_IN_TARGET_SITES, spelled out as 22 nt,
    or given as 'name=<22 nt>'.
    """
    if token in BUILT_IN_TARGET_SITES:
        return BUILT_IN_TARGET_SITES[token]

    name, sep, seq = token.partition("=")
    if not sep:
        name, seq = "", token

    seq = clean_sequence(seq, "target site")
    if len(seq) == 22:
        return name.strip(), seq

    known = ", ".join(BUILT_IN_TARGET_SITES)
    sys.exit(
        f"Error: '{token}' is neither a known target site ({known}) "
        f"nor a 22-nt DNA sequence."
    )


def _available(index: dict, kind: str) -> str:
    if not index:
        return f"no {kind} results are cached for this run"
    keys = sorted(index, key=lambda k: tuple(map(int, str(k).split("."))))
    return f"available {kind} sites: {', '.join(str(k) for k in keys)}"


def _parse_ami(token: str, rest: str, position: int, amirna_index: dict) -> dict:
    """'ami:N' from the cache, or 'ami:[name:]<21 nt>[:<21 nt star>]' typed out."""
    parts = [p.strip() for p in rest.split(":")] if rest.strip() else []

    if len(parts) == 1 and parts[0].isdigit():
        n = int(parts[0])
        if n not in amirna_index:
            sys.exit(f"Error: amiRNA result {n} does not exist — {_available(amirna_index, 'amiRNA')}.")
        entry = amirna_index[n]
        return amirna_module(f"amiRNA {n}", entry["amiRNA"], entry["amiRNA*"])

    if len(parts) == 1:
        return amirna_module(
            f"amiRNA {position}", clean_sequence(parts[0], "amiRNA"))

    if len(parts) == 2:
        return amirna_module(parts[0], clean_sequence(parts[1], "amiRNA"))

    if len(parts) == 3:
        return amirna_module(
            parts[0],
            clean_sequence(parts[1], "amiRNA"),
            clean_sequence(parts[2], "amiRNA*"))

    sys.exit(
        f"Error: '{token}' must look like 'ami:N' for a cached result, or "
        f"'ami:<name>:<21-nt amiRNA>[:<21-nt amiRNA*>]' for one typed out."
    )


def _parse_syn(token: str, rest: str, syn_index: dict) -> dict:
    """'syn:<target site>:<guides>', each guide cached or typed out."""
    site_token, _, refs = rest.partition(":")
    if not site_token.strip() or not refs.strip():
        sys.exit(
            f"Error: '{token}' must look like 'syn:<target site>:<sites>', "
            f"e.g. 'syn:NbmiR482a:1.1+2.1' or 'syn:NbmiR482a:NbDXS1=TAAACC...'."
        )
    target_name, target_site = _resolve_target_site(site_token.strip())

    guides = []
    for i, ref in enumerate(refs.split("+"), start=1):
        ref = ref.strip()
        if not ref:
            continue

        name, sep, seq = ref.partition("=")
        if sep:
            guides.append((name.strip(), clean_sequence(seq, f"syn-tasiRNA {i}")))
            continue

        if re.fullmatch(r"\d+\.\d+", ref):
            if ref not in syn_index:
                sys.exit(
                    f"Error: syn-tasiRNA site '{ref}' does not exist — "
                    f"{_available(syn_index, 'syn-tasiRNA')}."
                )
            guides.append((f"syn-tasiRNA {ref}", syn_index[ref]))
            continue

        guides.append((f"syn-tasiRNA {i}", clean_sequence(ref, f"syn-tasiRNA {i}")))

    return syntasirna_module(target_name, target_site, guides)


def parse_modules(spec: str, amirna_index: dict = None, syn_index: dict = None) -> list:
    """
    Turn an ordered, comma-separated module specification into the module
    list build_insert() expects.

    Every sequence is either addressed into a previous run's optimal results,
    the way the rest of P-SAMS addresses them, or typed out in full so the
    tool can be used on its own with no run behind it.

    Tokens, in the order the modules are to be laid down:

      ami:N                           cached optimal amiRNA result N
      ami:<name>:<ami>[:<star>]       an amiRNA typed out; the amiRNA* is
                                      derived from it when left out
      mir173a[:<name>]                the fixed AtMIR173a precursor
      syn:<site>:<guides>             a syn-tasiRNA module, where <site> is a
                                      built-in target site name, 22 nt, or
                                      'name=<22 nt>', and <guides> is a
                                      '+'-joined list whose items are each a
                                      cached 'geneset.site' reference, a
                                      'name=<21 nt>' pair, or bare 21 nt

    For example: 'ami:1,syn:NbmiR482a:1.1+2.1'
             or: 'ami:amiR-NbSu:TGTATGACTCCCGGAATTCCA,mir173a'
    """
    amirna_index = amirna_index or {}
    syn_index = syn_index or {}
    modules = []

    for raw in spec.split(","):
        token = raw.strip()
        if not token:
            continue

        kind, _, rest = token.partition(":")
        kind = kind.strip().lower()

        if kind in ("mir173a", "atmir173a"):
            modules.append(mir173a_module(rest.strip() or "AtMIR173a"))
        elif kind == "ami":
            modules.append(_parse_ami(token, rest, len(modules) + 1, amirna_index))
        elif kind == "syn":
            modules.append(_parse_syn(token, rest, syn_index))
        else:
            sys.exit(
                f"Error: '{token}' is not a module. Expected 'ami:...', "
                f"'syn:...' or 'mir173a'."
            )

    if not modules:
        sys.exit("Error: --modules did not name a single module.")
    return modules


def spec_uses_cache(spec: str) -> bool:
    """
    Whether a module specification refers to a previous run's results at all.
    A fully typed-out specification needs no cached run, and no output folder
    to find one in.
    """
    for raw in spec.split(","):
        token = raw.strip()
        kind, _, rest = token.partition(":")
        kind = kind.strip().lower()
        if kind == "ami" and rest.strip().isdigit():
            return True
        if kind == "syn":
            _, _, refs = rest.partition(":")
            for ref in refs.split("+"):
                if "=" not in ref and re.fullmatch(r"\d+\.\d+", ref.strip()):
                    return True
    return False


# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

def safe_name(name: str, fallback: str = "construct") -> str:
    """Make a construct name usable as a filename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name or fallback)


def design_basename(design: dict) -> str:
    """The stem the HTML gives this design's downloads."""
    return f"{safe_name(design['construct_name'], design['name_fallback'])}_{design['file_tag']}"


def _module_block(module: dict, i: int, shape: str) -> str:
    """One module's block in the TXT report, in the shape the HTML uses."""
    if shape == "ami":
        return "\n".join([
            f"  Module {i}: {module['name']}",
            f"    amiRNA  (21nt): {module['amiRNA']}",
            f"    amiRNA* (21nt): {module['amiRNA*']}",
        ])

    if shape == "syn":
        lines = [f"  Module {i}"]
        lines.append(f"    Target site \"{module['target_name'] or '-'}\" (22nt): {module['target_site']}")
        for j, (name, seq) in enumerate(module["guides"], start=1):
            lines.append(f"    syn-tasiRNA {j} \"{name}\" (21nt): {seq}")
        return "\n".join(lines)

    if module["type"] == "AtMIR173a":
        return "\n".join([
            f"  Module {i} (AtMIR173a foldback precursor)",
            f"    \"{module['name']}\" ({len(module['sequence'])}nt): {module['sequence']}",
        ])
    if module["type"] == "amiRNA":
        return "\n".join([
            f"  Module {i} (amiRNA)",
            f"    amiRNA \"{module['name']}\" (21nt): {module['amiRNA']}",
            f"    amiRNA* (21nt): {module['amiRNA*']}",
        ])
    lines = [f"  Module {i} (syn-tasiRNA)"]
    lines.append(f"    Target site \"{module['target_name'] or '-'}\" (22nt): {module['target_site']}")
    for j, (name, seq) in enumerate(module["guides"], start=1):
        lines.append(f"    syn-tasiRNA {j} \"{name}\" (21nt): {seq}")
    return "\n".join(lines)


def build_txt(design: dict) -> str:
    shape = "bc" if design["vector"].startswith("pMDC32B") else design_shape(design["modules"])
    blocks = [
        _module_block(m, i, shape)
        for i, m in enumerate(design["modules"], start=1)
    ]
    return "\n".join([
        f"Construct name: {design['construct_name']}",
        f"Type: {design['type']}",
        "",
        "Modules:",
        "\n\n".join(blocks),
        "",
        f"Full {design['strandedness']} insert (5'->3'): {design['insert']}",
        "",
        "Note: Verify against your exact construct design before ordering.",
    ])


def build_fasta(design: dict) -> str:
    name = safe_name(design["construct_name"], design["name_fallback"])
    return f">{name}_full_insert\n{design['insert']}"


def build_csv(design: dict) -> str:
    shape = "bc" if design["vector"].startswith("pMDC32B") else design_shape(design["modules"])
    rows = [
        ("construct_name", design["construct_name"]),
        ("type", design["csv_type"]),
    ]

    for i, module in enumerate(design["modules"], start=1):
        if shape == "ami":
            rows.append((f"module_{i}_name", module["name"]))
            rows.append((f"module_{i}_amiRNA", module["amiRNA"]))
            rows.append((f"module_{i}_amiRNA_star", module["amiRNA*"]))
            continue
        if shape == "syn":
            rows.append((f"module_{i}_target_name", module["target_name"]))
            rows.append((f"module_{i}_target_seq", module["target_site"]))
            for j, (name, seq) in enumerate(module["guides"], start=1):
                rows.append((f"module_{i}_syn_{j}_name", name))
                rows.append((f"module_{i}_syn_{j}_seq", seq))
            continue

        if module["type"] == "AtMIR173a":
            rows.append((f"module_{i}_type", "AtMIR173a_precursor"))
            rows.append((f"module_{i}_name", module["name"]))
            rows.append((f"module_{i}_sequence", module["sequence"]))
        elif module["type"] == "amiRNA":
            rows.append((f"module_{i}_type", "amiRNA"))
            rows.append((f"module_{i}_amiRNA_name", module["name"]))
            rows.append((f"module_{i}_amiRNA", module["amiRNA"]))
            rows.append((f"module_{i}_amiRNA_star", module["amiRNA*"]))
        else:
            rows.append((f"module_{i}_type", "syn-tasiRNA"))
            rows.append((f"module_{i}_target_name", module["target_name"]))
            rows.append((f"module_{i}_target_seq", module["target_site"]))
            for j, (name, seq) in enumerate(module["guides"], start=1):
                rows.append((f"module_{i}_syn_{j}_name", name))
                rows.append((f"module_{i}_syn_{j}_seq", seq))

    rows.append(("full_insert_5to3", design["insert"]))

    buf = io.StringIO()
    writer = csv.writer(buf, quoting=csv.QUOTE_ALL, lineterminator="\n")
    writer.writerow(["field", "value"])
    writer.writerows(rows)
    return buf.getvalue().rstrip("\n")


def write_exports(design: dict, directory) -> list:
    """
    Write this design's TXT, FASTA and CSV companions into `directory`, named
    the way the Oligo Designer Suite names its downloads, and return what was
    written.
    """
    base = design_basename(design)
    written = []
    for filename, text in (
        (f"{base}_design.txt", build_txt(design)),
        (f"{base}_design.fasta", build_fasta(design)),
        (f"{base}_summary.csv", build_csv(design)),
    ):
        path = directory / filename
        path.write_text(text + "\n")
        written.append(path)
    return written
