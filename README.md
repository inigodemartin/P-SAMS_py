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

For `syntasiRNA`, you can define multiple groups of genes/sequences by
separating the groups with `:`, e.g. `-a gen1,gen2:gen3,gen4`.

### Generated results

```
runs/Nbe01g01610/Nbe01g01610.7_psams_output/
├── Nbe01g01610.7_optimal_results.tsv      # optimal amiRNAs (no off-targets)
├── Nbe01g01610.7_suboptimal_results.tsv   # suboptimal amiRNAs (with off-targets)
├── Nbe01g01610.7_psams.json               # final combined result
└── tf_results/
    ├── site_0001_TargetFinder_result.json # TargetFinder output per candidate
    └── ...
```

If `-o` already contains complete results for those accessions, the script
detects this and exits without re-running the analysis.

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
