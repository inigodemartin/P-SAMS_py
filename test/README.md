# Tests

Regression tests for the bugs reported by Alberto after running P-SAMS
against the pepper (*Capsicum annuum*, 34,899 transcripts) transcriptome.
They use only the Python standard library plus `perl`, and build their own
synthetic FASTA and SQLite inputs in a temporary directory, so no external
database or downloaded transcriptome is needed. The whole suite runs in a
couple of seconds.

## Running

From the repository root:

```bash
conda activate p-sams-py
python3 test/test_pepper_bugs.py
```

Expected output: `Ran 4 tests` and `OK`.

## What each test covers

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
