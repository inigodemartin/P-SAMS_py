import sys
import subprocess
import re

def convert_fasta_to_string(fasta_file):
    fasta_str = ''
    with open(fasta_file, 'r') as f:
        for line in f:
            fasta_str += line
    return fasta_str
        

def build_fg_index_fasta(fasta: str, debug: bool = False) -> dict:
    """
    Parse FASTA-formatted string and build a dictionary:
    {sequence_id: sequence}

    Parameters
    ----------
    fasta : str
        FASTA formatted text
    debug : bool
        Print debug information

    Returns
    -------
    dict
        Dictionary mapping IDs to sequences
    """

    ids = {}
    current_id = None

    if debug:
        print("Building foreground index...", end="", file=sys.stderr)

    for line in fasta.splitlines():
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Header line
        if line.startswith(">"):
            current_id = line[1:]
            ids[current_id] = ""

        # Sequence line
        else:
            if current_id is None:
                raise ValueError("FASTA format error: sequence without header")
            ids[current_id] += line

    if debug:
        print(f"{len(ids)} sequences loaded...done", file=sys.stderr)

    return ids


def build_fg_index(accessions, conn, species, mRNA_fa, min_length=9, debug=False):
    """
    Build foreground index from gene accessions using SQLite + samtools.

    Parameters
    ----------
    accessions : list[str]
        List of gene/transcript IDs
    conn : sqlite3.Connection
        Active database connection
    species : str
        Species name (used for ID normalization)
    mRNA_fa : str
        Path to mRNA FASTA
    min_length : int
        Minimum length of accession
    debug : bool
        Debug flag

    Returns
    -------
    dict
        {transcript_id: sequence}
    """

    import re
    import subprocess
    import sys
    import sqlite3
    ids = {}
    cursor = conn.cursor()

    if debug:
        print("Building foreground index...", end="", file=sys.stderr)

    # SQLite uses ? instead of %s
    query = "SELECT * FROM annotation WHERE transcript LIKE ?"

    for accession in accessions:

        # --------------------------------------------------
        # Normalize accession (UNCHANGED LOGIC)
        # --------------------------------------------------
        accession = re.sub(r"\.\d{1,2}$", "", accession)

        if species == "Z_MAYS":
            if accession.startswith("AC"):
                accession = accession.replace("FGT", "FG")
            else:
                accession = re.sub(r"_T\d{2}$", "", accession)

        elif species == "C_REINHARDTII":
            accession = re.sub(r"\.t\d$", "", accession)

        # Length check
        if len(accession) < min_length:
            print(f"Gene {accession} not found in database!", file=sys.stderr)
            sys.exit(1)

        # --------------------------------------------------
        # Query SQLite database
        # --------------------------------------------------
        exists = False
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")

        cursor.execute(query, (f"{accession}%",))
        results = cursor.fetchall()

        for result in results:
            exists = True
            transcript = result[0] if isinstance(result, tuple) else result["transcript"]
            ids[transcript] = ""

            # --------------------------------------------------
            # Run samtools faidx (UNCHANGED)
            # --------------------------------------------------
            cmd = ["samtools", "faidx", mRNA_fa, transcript]
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                text=True
            )

            for line in process.stdout:
                if line.startswith(">"):
                    continue
                ids[transcript] += line.strip()

            process.stdout.close()
            process.wait()

        if not exists:
            print(f"Gene {accession} not found in database!", file=sys.stderr)
            sys.exit(1)

    if debug:
        print("done", file=sys.stderr)

    return ids

def parse_list(flat_list: str, sep: str) -> list:
    """
    Parse a delimited string into a list.

    Parameters
    ----------
    flat_list : str
        Input string (e.g. "gene1,gene2,gene3")
    sep : str
        Delimiter character

    Returns
    -------
    list
        List of elements
    """

    if sep in flat_list:
        return flat_list.split(sep)
    else:
        return [flat_list]