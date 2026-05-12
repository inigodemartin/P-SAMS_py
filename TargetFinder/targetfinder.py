#!/usr/bin/env python3

import os
import re
import sys
import tempfile
import subprocess
import argparse

DEBUG = False

################################################################################
# Variable declarations
################################################################################

TMPDIR = os.environ.get("TMPDIR", "/tmp")
FASTA_BIN = "ssearch35_t"

if DEBUG:
    LOG = open("targetfinder.log", "w")

def log(msg):
    if DEBUG:
        LOG.write(msg + "\n")

################################################################################
# Argument parsing (equivalent to getopts)
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("-s", required=True)
    parser.add_argument("-d", required=True)
    parser.add_argument("-q", default="query")
    parser.add_argument("-c", type=float, default=4)
    parser.add_argument("-t", type=int, default=1)
    parser.add_argument("-p", default="classic")
    parser.add_argument("-r", action="store_true")
    parser.add_argument("-h", action="store_true")

    args = parser.parse_args()

    if args.h:
        var_error()

    if args.p not in ["classic", "gff", "json", "table"]:
        print(f"Output format {args.p} is not valid.", file=sys.stderr)
        sys.exit(1)

    return args

################################################################################
# FASTA execution
################################################################################

def fasta(input_file, db, threads):
    cmd = [
        FASTA_BIN,
        "-n", "-H", "-Q",
        "-f", "-16",
        "-r", "+15/-10",
        "-g", "-10",
        "-w", "100",
        "-W", "25",
        "-E", "100000",
        "-i", "-U",
        "-T", str(threads),
        input_file,
        db,
        "1"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.splitlines()

################################################################################
# FASTA parser
################################################################################

def fasta_parser(lines, name):
    output = []

    step = 0
    hit_accession = ""
    miRNA = ""
    target = ""
    spacer_count = 0
    length = 0
    short_id = ""

    for line in lines:
        if step == 0:
            m = re.match(r"^>>(.+)\(\d+ nt\)", line)
            if m:
                hit_accession = m.group(1).rstrip()
                short_id = hit_accession[:6]
                step = 1

        elif step == 1:
            m = re.search(rf"({name}-\s+)(\D+)", line)
            if m:
                spacer = m.group(1)
                miRNA = m.group(2).replace(" ", "")
                spacer_count = len(re.findall(r"[A-Za-z0-9\-_ ]", spacer))
                length = len(miRNA)
                miRNA = miRNA.translate(str.maketrans("AUGC", "UACG"))

            elif short_id in line:
                m2 = re.search(rf".{{{spacer_count}}}(\D{{{length}}})", line)
                if m2:
                    target = m2.group(1).replace(" ", "")
                    output.append((hit_accession, miRNA, target))
                    step = 0

    return output

################################################################################
# Base pairing scoring
################################################################################

bp = {
    "AU":0,"UA":0,"GC":0,"CG":0,
    "GU":0.5,"UG":0.5,
    "AC":1,"CA":1,"AG":1,"GA":1,
    "UC":1,"CU":1,
    "A-":1,"U-":1,"G-":1,"C-":1,
    "-A":1,"-U":1,"-G":1,"-C":1,
    "AA":1,"UU":1,"CC":1,"GG":1
}

def bp_score(parsed, cutoff, query_name):
    results = []

    for hit_accession, miR_seq, target_seq in parsed:

        if not re.match(r"^[AGCU-]+$", target_seq):
            continue
        if len(miR_seq) != len(target_seq):
            continue
        if re.search(r"-.*-", miR_seq):
            continue
        if re.search(r"-.*-", target_seq):
            continue
        if "-" in miR_seq and "-" in target_seq:
            continue

        mismatch = 0
        gu = 0
        score = 0
        homology = ""

        miRNA = list(miR_seq)
        target = list(target_seq)

        for cycle in range(1, len(miRNA)+1):
            miR_base = miRNA.pop()
            target_base = target.pop()

            val = bp.get(miR_base + target_base, 1)

            if cycle > 1 and cycle <= 13:
                val *= 2

            if val == 1 or val == 2:
                mismatch += 1
                homology += " "
            elif val == 0.5 or val == 1:
                gu += 1
                homology += "."
            else:
                homology += ":"

            score += val

        total = mismatch + gu

        if total > 7 or mismatch > 4 or gu > 4:
            continue

        if score <= cutoff:
            results.append({
                "hit_accession": hit_accession,
                "miR_seq": miR_seq,
                "query_name": query_name,
                "target_seq": target_seq,
                "score": score,
                "homology_string": homology[::-1],
                "target_start": 0,
                "target_end": 0,
                "target_strand": 0
            })

    return results

################################################################################
# Output formats
################################################################################

def print_classic(targets):
    for t in targets:
        print(f"query={t['query_name']}, target={t['hit_accession']}, score={t['score']}, "
              f"range={t['target_start']}-{t['target_end']}, strand={t['target_strand']}\n")
        print(f"target  5' {t['target_seq']} 3'")
        print(f"           {t['homology_string']}")
        print(f"query   3' {t['miR_seq']} 5'\n")

def print_table(targets):
    for t in targets:
        strand = "+" if t["target_strand"] == 1 else "-"
        print("\t".join(map(str, [
            t["query_name"],
            t["hit_accession"],
            t["target_start"],
            t["target_end"],
            strand,
            t["score"],
            t["target_seq"],
            t["homology_string"],
            t["miR_seq"]
        ])))

################################################################################
# Error help
################################################################################

def var_error():
    print("""
TargetFinder: Plant small RNA target prediction tool.

Usage: targetfinder.py -s <sequence> -d <database>
""", file=sys.stderr)
    sys.exit(1)

################################################################################
# MAIN
################################################################################

def main():
    args = parse_args()

    sRNA = args.s.upper().replace("T", "U")
    query_name = args.q

    name = query_name[:5]

    with tempfile.NamedTemporaryFile(delete=False, dir=TMPDIR, mode="w") as f:
        f.write(f">{query_name}\n{sRNA}")
        infile = f.name

    fasta_out = fasta(infile, args.d, args.t)
    parsed = fasta_parser(fasta_out, name)
    targets = bp_score(parsed, args.c, query_name)

    if targets:
        targets.sort(key=lambda x: x["score"])

        if args.p == "classic":
            print_classic(targets)
        elif args.p == "table":
            print_table(targets)

    else:
        print(f"No results for {query_name}")

    os.unlink(infile)

if __name__ == "__main__":
    main()