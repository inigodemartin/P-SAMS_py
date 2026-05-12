#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO


def format_fasta(in_fasta, species, version, output_folder):
    """
    Parse FASTA using Bio.SeqIO:
    - Keep only transcript ID
    - Remove sequences with ChrSy.fgenesh / ChrUn.fgenesh
    """

    out_fasta = output_folder / f"{species}.{version}.transcripts.fasta"

    if out_fasta.exists() and out_fasta.stat().st_size > 0:
        return out_fasta

    with open(out_fasta, "w") as outfile:
        for record in SeqIO.parse(in_fasta, "fasta"):

            header = record.description

            # Filter unwanted entries
            if "ChrSy.fgenesh" in header or "ChrUn.fgenesh" in header:
                continue

            # Keep only first token of header
            record.id = record.id.split()[0]
            record.description = record.id

            SeqIO.write(record, outfile, "fasta")

    return out_fasta


def format_annot_info(annot_file, species, version, output_folder):
    """
    Parse annotation file and generate transcript-to-description table.
    """

    outfile = output_folder / f"{species}.{version}.annotation.txt"
    if outfile.exists() and outfile.stat().st_size > 0:
        return outfile

    with open(outfile, "w") as out, open(annot_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")

            if len(fields) < 16:
                continue

            (
                _id, gene, transcript, protein,
                pfam, panther, kog, ec, ko, go,
                best_ath, ath_sym, ath_desc,
                best_os, os_sym, os_desc
            ) = fields

            # -------------------------
            # Description
            # -------------------------
            description = ""

            if ath_desc:
                description += ath_desc
                if ath_sym:
                    description += f" ({ath_sym})"
            elif os_desc:
                description += os_desc
                if os_sym:
                    description += f" ({os_sym})"
            else:
                description = "unknown"

            # -------------------------
            # CONTAINS
            # -------------------------
            contains = []

            if any([pfam, panther, kog, ec, ko, go]):
                if pfam:
                    contains.append(f"PFAM:{pfam}")
                if panther:
                    contains.append(f"PANTHER:{panther}")
                if kog:
                    contains.append(f"KOG:{kog}")
                if ec:
                    contains.append(f"EC:{ec}")
                if ko:
                    contains.append(f"KEGG:{ko}")
                if go:
                    contains.append(f"GO:{go}")

                description += "; CONTAINS " + "; ".join(contains)

            # -------------------------
            # BEST
            # -------------------------
            best = []

            if best_ath:
                best.append(best_ath)
            if best_os:
                best.append(best_os)

            if best:
                description += "; BEST:" + ",".join(best)

            out.write(f"{transcript}\t{description}\n")

    return outfile


def parse_args():
    parser = argparse.ArgumentParser(
        prog="phytozome2psams.py",
        description="Parse Phytozome transcript FASTA and annotation files."
    )

    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-a", "--annotation", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-v", "--version", required=True)

    return parser.parse_args()



def phytozome2psams(fasta, species, version, annotation, output_folder):
    """
    Function for running the script as a function
    """
    out_fasta = format_fasta(fasta, species, version, output_folder)
    out_annot = format_annot_info(annotation, species, version, output_folder)
    return out_fasta, out_annot

