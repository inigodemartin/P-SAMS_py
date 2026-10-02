# Tests

Two suites, both runnable on their own. Neither needs a transcriptome, an
annotation database or any network access.

| File | Covers |
|---|---|
| `test_pepper_bugs.py` | the four bugs Alberto reported from the pepper transcriptome |
| `test_insert_design.py` | the full-insert vectors PVX, TRV and the B/c multimodule entry |

```bash
conda activate p-sams-py
python3 test/test_pepper_bugs.py
python3 test/test_insert_design.py
```

## `test_pepper_bugs.py`

Regression tests for the bugs reported by Alberto after running P-SAMS
against the pepper (*Capsicum annuum*, 34,899 transcripts) transcriptome.
They use only the Python standard library plus `perl`, and build their own
synthetic FASTA and SQLite inputs in a temporary directory, so no external
database or downloaded transcriptome is needed. The whole suite runs in a
couple of seconds.

Expected output: `Ran 9 tests` and `OK`.


The resume tests are pure unit tests over `load_resume_state` and need no
subprocess at all.

| Test | Bug it guards against |
|---|---|
| `test_database_path_with_spaces_and_parentheses` | `TargetFinder/targetfinder.pl`, `sub fasta()` — an unquoted database path broke the `ssearch36` shell call for paths holding spaces or parentheses (a synced Google Drive folder, say), and the pipeline reported a false "0 optimal sites" |
| `test_targetfinder_json_survives_backslash_in_header` | `TargetFinder/targetfinder.pl`, `sub print_json()` — FASTA headers containing a backslash (two pepper genes do) produced invalid JSON, since `\+` is not a valid JSON escape |
| `test_off_target_check_json_survives_backslash_in_description` | `src/utils.py`, `off_target_check()` — the same backslash problem by a second, independent route, through the annotation description. Only reachable with `-u` |
| `test_fasta_mode_suboptimal_candidate_writes_valid_json` | `src/pipeline.py`, `serial_jobs()` — the `-f` plus suboptimal-candidate branch stored the TargetFinder result as a list rather than joined text, raising a `TypeError` whenever off-targets were checked |

## Synthetic data

There is no committed data file. Every test builds what it needs from two
constants at the top of the module: a 21-nt guide and a transcript carrying
its exact reverse complement, which guarantees TargetFinder a perfect hit
without real sequence data. The suboptimal-branch test writes that same
transcript twice, under two accessions, so the second one is an unavoidable
off-target and the candidate is forced down the suboptimal branch.

## `test_insert_design.py`

Covers the vectors whose output is a whole synthetic fragment rather than a
pair of annealing oligos, ported from Alberto's Oligo Designer Suite HTML.
Every reference value is the HTML's own: its literal flank and scaffold
sequences, the module and insert lengths its footer states, and its worked
example `amiR-NbSu`.

What it pins down:

- The PVX, TRV and B/c multimodule flanks, base for base.
- One BsaI site per strand in the B/c multimodule insert, and none at all in
  the viral ones.
- An amiRNA module of 89 nt, and the insert lengths the HTML reports for
  every one of its worked examples.
- The B/c multimodule example reproduced byte for byte.
- That P-SAMS's own amiRNA\* convention builds a byte-identical module to
  the differently-split convention the HTML's TRV app uses. The two label
  the star differently but must never produce different DNA.
- That modules are butt-joined, with nothing inserted between them.
- Module specifications resolving against a run's cached results, including
  a hybrid, and every malformed or dangling reference being refused.

Standalone use, with no cached run behind the design:

- The amiRNA\* derived from the amiRNA, against both the HTML's worked
  example and its TRV app's independent derivation.
- Sequence cleaning and the length checks, matching the HTML's own.
- Cached and typed-out guides mixing inside one module, named target sites,
  and a renamed AtMIR173a module.
- `clone_vector.py` building an insert in a directory holding no run at all.

HTML parity:

- The `Type:` label for every vector and module mix.
- An AtMIR173a-only insert no longer being described as a MIR390 amiRNA.
- The B/c insert reported as single-stranded and the viral ones as
  double-stranded.
- Export filenames following the HTML's downloads, fallbacks included.

And the foldback guard: a run cached as `monocot` is refused rather than
wrapped in the eudicot scaffold, both as a unit check and through
`clone_vector.py`.

Expected output: `Ran 39 tests` and `OK`.
