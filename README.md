# P-SAMS_py

Python version of **P-SAMS** (Plant Small RNA Maker Suite): a tool for
designing **amiRNA** (artificial microRNA) and **syn-tasiRNA**, also predicting
possible targets and off-targets via **TargetFinder**.

## Requirements

- [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
- Python 3.9+ (tested with 3.12), [Biopython](https://biopython.org/), Perl 5
  and `samtools` are all installed through the conda environment (see below)

The `ssearch36` binaries (Linux, macOS arm64/x86_64) are already included in
`TargetFinder/` and should have execute permissions. If not:

```bash
chmod +x TargetFinder/ssearch36_* TargetFinder/targetfinder.pl
```

## Installation

```bash
git clone https://github.com/inigodemartin/P-SAMS_py.git
cd P-SAMS_py
conda env create -f environment.yml
conda activate p-sams-py
```

No additional "build" step is required: `psams.py` and `create_db.py` are run
directly with `python3` once the environment is activated.

To update an already-created environment after `environment.yml` changes:

```bash
conda env update -f environment.yml --prune
```

---

## 1. Building a species database

Before you can analyze genes of a species with `-a/--accessions`, you need to
generate its k-mer database (used to predict off-targets) and register it in
`psams.conf`.

### Input files

You need the typical **Phytozome** files for your species:

- A transcript (mRNA) FASTA, e.g. `Nicotiana_benthamiana.transcript.fa`
- An annotation file `*.annotation_info.txt` (Phytozome format, 16
  tab-separated columns)

### Command

```bash
python3 create_db.py \
    -f Nicotiana_benthamiana.transcript.fa \
    -d Nicotiana_benthamiana.annotation_info.txt \
    -s Nicotiana_benthamiana \
    -v v1.0 \
    -k 15
```

Arguments:

| Flag | Description |
|------|-------------|
| `-f, --fasta`        | Transcript FASTA (Phytozome format) |
| `-d, --descriptions` | Annotation file (`annotation_info.txt`) |
| `-s, --species`      | Species name/code (used later with `-s` in `psams.py`) |
| `-v, --version`      | Transcriptome version (only used to name the generated files) |
| `-k, --ksize`        | K-mer size for the off-target database. **Must be 15** (matches the `SEED` used by `psams.py`) |
| `-l, --lowmem`       | Optional. Builds the database directly, slower but using much less RAM |

This step can take a while (it iterates over the whole transcriptome generating k-mers).

### What it generates

- `db/<species>_<version>_prueba.db`: SQLite database with the `kmers` and
  `annotation` tables.
- Intermediate formatted files alongside the input FASTA
  (`<species>.<version>.transcripts.fasta`, `<species>.<version>.annotation.txt`)
  and its `.fai` index.
- A new entry in `psams.conf` (in the project root) with the format:

```ini
[Nicotiana_benthamiana]
mRNA=/absolute/path/to/Nicotiana_benthamiana.transcript.fa
sql=/absolute/path/to/project/db/Nicotiana_benthamiana_v1.0_prueba.db
```

`psams.conf` is not version-controlled (it's in `.gitignore`); each
installation has its own with the species it has built.

---

## 2. Running an analysis

Once the species is registered in `psams.conf`, you can analyze a gene by its
accession:

```bash
python3 psams.py -a Nbe01g01610.7 -s Nicotiana_benthamiana -o runs/Nbe01g01610 -u
```

### Main flags

| Flag | Description |
|------|-------------|
| `-a, --accessions`   | Gene accession(s), comma-separated. Requires `-s` |
| `-f, --fasta`        | Alternative to `-a`: FASTA file with the sequence(s) to analyze |
| `-s, --species`      | Species as it appears in `psams.conf`. Required if using `-a`, or if using `-f` and off-target prediction is wanted |
| `-o, --output_path`  | Folder where results are created (defaults to the current directory) |
| `-c, --construct`    | `amiRNA` (default) or `syntasiRNA` |
| `-t, --foldback`     | `eudicot` (default) or `monocot` |
| `-n, --noofftarget`  | Disables off-target prediction with TargetFinder |
| `-u, --unlimit`      | Don't limit to 3 optimal results: go through all possible candidates (slower) |
| `-j, --jobs`         | Number of TargetFinder jobs to run in parallel. Default = 1 (serial) |
| `-V, --vector`       | Prompts an interactive menu to pick a cloning vector (see below). Adds vector-specific cloning oligos to the output |
| `-T, --target-site`  | 22-nt miRNA target site sequence. Only required when the `pMDC32B-B/c` vector is selected. If omitted and needed, it is requested interactively |

### Cloning vectors (`-V`)

When `-V` is passed, the script prompts an interactive menu listing the vectors
available for the chosen `--construct`, so you just pick a number instead of
typing the exact vector name. The cloning oligos in the output JSON are then
computed using the selected vector's architecture. The available vectors
depend on the construct type:

**amiRNA** (`-c amiRNA`):

| Vector | Overhang |
|--------|----------|
| `pMDC32B-AtMIR390a-B/c` | TGTA / AATG |
| `pMDC32B-OsMIR390-B/c` | CTTG / CATG |
| `pMDC32B-BS-AtMIR390a-B/c` | TGTA / AATG |
| `pMDC32B-BS-AtMIR390a-A18G-B/c` | TGTG / AATG |

**syn-tasiRNA** (`-c syntasiRNA`):

| Vector | Overhang | Notes |
|--------|----------|-------|
| `pMDC32B-B/c` | TGTA / AATG | Requires `-T` (22-nt miRNA target site) |
| `pMDC32B-AtTAS1c-B/c` | ATTA / GTTC | |
| `pMDC32B-AtTAS1c-D2-B/c` | TTTA / CCGA | |
| `pMDC32B-AtmiR173aTS-B/c` | TTTA / CCGA | |
| `pMDC32B-NbmiR482aTS-B/c` | TTTA / CCGA | |
| `pMDC32B-SlmiR482bTS-B/c` | TTTA / CCGA | |

**Example — amiRNA with a cloning vector:**

```bash
python3 psams.py -a Nbe01g01610.7 -s Nicotiana_benthamiana \
    -o runs/Nbe01g01610 -V
```

```
Available cloning vectors:
  1. pMDC32B-AtMIR390a-B/c
  2. pMDC32B-OsMIR390-B/c
  3. pMDC32B-BS-AtMIR390a-B/c
  4. pMDC32B-BS-AtMIR390a-A18G-B/c
Select a vector [1-4]: 1
```

The output filename gets the chosen vector appended (e.g.
`Nbe01g01610.7_pMDC32B-AtMIR390a-Bc_psams.json`), and the JSON itself
includes a top-level `"vector"` key plus vector-specific
`"Forward Oligo"` / `"Reverse Oligo"` per candidate.

**Example — syn-tasiRNA with target site:**

```bash
python3 psams.py -a Nbe01g01610.7 -s Nicotiana_benthamiana \
    -c syntasiRNA -o runs/Nbe01g01610 -V
```

If the vector you pick from the menu is `pMDC32B-B/c` and `-T` wasn't
given on the command line, you'll be prompted for the 22-nt target site
sequence interactively. The `_psams.json` output filename and content
follow the same convention: vector name in the filename, plus a
`"cloning_oligos"` key with the final forward and reverse oligos ready
to order.

For `syntasiRNA`, you can define multiple groups of genes/sequences by
separating the groups with `:`, e.g. `-a gen1,gen2:gen3,gen4`.

**Re-running with a different vector:** every run writes a vector-agnostic
`{accession}_psams.json` cache alongside the vector-specific file. If you
run the same input again with `-V` and a full run's results already exist,
P-SAMS reuses that cache instead of re-running TargetFinder, and only
(re)computes the cloning oligos for the newly selected vector:

```
P-SAMS already executed for Nbe01g01610.7.
Reusing cached results and generating cloning oligos for vector 'pMDC32B-OsMIR390-B/c'.
Output: runs/Nbe01g01610/Nbe01g01610.7_psams_output/Nbe01g01610.7_pMDC32B-OsMIR390-Bc_psams.json
Exiting script.
```

### Generated results

```
runs/Nbe01g01610/Nbe01g01610.7_psams_output/
├── Nbe01g01610.7_optimal_results.tsv               # optimal amiRNAs (no off-targets)
├── Nbe01g01610.7_suboptimal_results.tsv            # suboptimal amiRNAs (with off-targets)
├── Nbe01g01610.7_psams.json                        # vector-agnostic cache / final result
├── Nbe01g01610.7_pMDC32B-OsMIR390-Bc_psams.json    # only created when -V is used
└── tf_results/
    ├── site_0001_TargetFinder_result.json # TargetFinder output per candidate
    └── ...
```

If `-o` already contains complete results for those accessions, the script
detects this and exits without re-running the analysis — unless `-V` is
used, in which case it reuses the cached `_psams.json` to generate the
vector-specific file without redoing the (slow) TargetFinder search.

---

## 3. Quick test (without a database)

To check that the installation works without building any database, you can
pass your own sequence with `-f` and disable off-target prediction with `-n`
(this mode doesn't use Perl, samtools, or `psams.conf`).

Create a file `example.fasta`:

```fasta
>my_transcript
ATGGCGGATTCAGAGAAGCCGGTTACCGGAAGCTTGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTG
```

And run it:

```bash
python3 psams.py -f example.fasta -n -o runs/quick_test
```

This will generate `runs/quick_test/my_transcript_psams_output/` with up to 3
optimal amiRNA designs (guide, star, cloning oligos) computed directly from the
input sequence, without checking off-targets.
