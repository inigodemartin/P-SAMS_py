#!/usr/bin/env python3
import sys
import sqlite3

site_length = 21

def connect_sqlite(db_file):
    """Connect to SQLite database."""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("Connected to SQLite database")
    return conn, cursor


def create_tables(cursor):
    """Create required tables."""
    
    cursor.execute("""
        DROP TABLE IF EXISTS kmers;
    """)
    
    cursor.execute("""
        CREATE TABLE kmers (
            kmer TEXT PRIMARY KEY,
            transcripts TEXT NOT NULL
        );
    """)

    cursor.execute("""
        DROP TABLE IF EXISTS annotation;
    """)

    cursor.execute("""
        CREATE TABLE annotation (
            transcript TEXT PRIMARY KEY,
            description TEXT
        );
    """)

    print("Tables created")


def process_fasta_load(fasta, ksize, conn):
    offset = site_length - ksize - 1
    transcript = None
    seq = ""
    count = 0
    cursor = conn.cursor()

    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if transcript is not None:
                    kmerize_load(seq,transcript,cursor,ksize,offset)
                    count += 1
                    if count % 100 == 0:
                        conn.commit()
                    print("Processed transcripts for kmer database: "+ str(count),end="\r",file=sys.stderr)
                transcript = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        # Process last transcript
        if transcript is not None:
            kmerize_load(seq,transcript,cursor,ksize,offset)
            count += 1
    conn.commit()
    print("\nProcessed " + str(count),file=sys.stderr)


def kmerize_load(seq, transcript, cursor, ksize, offset):

    seq = seq.upper()
    length = len(seq)
    seen = set()

    for i in range(length - site_length + 1):
        window = seq[i:i + site_length]
        kmer = window[offset:offset + ksize]
        
        # Avoid processing same kmer if it is more than once in the same transcript
        if kmer in seen:
            continue

        seen.add(kmer)
        cursor.execute(
            """
            SELECT transcripts
            FROM kmers
            WHERE kmer = ?
            """,
            (kmer,))
        row = cursor.fetchone()
        # Kmer does not exist yet
        if row is None:
            cursor.execute(
                """
                INSERT INTO kmers
                (kmer, transcripts)
                VALUES (?, ?)
                """,
                (kmer, transcript)
            )

        # Kmer already exists
        else:
            # se te ocurre alguna forma de hacerlo mas rapido?? porque cuando entra en este else se ralentiza bastante
            # y claro a medida de que la base de datos se va llenando, cada vez va entrando mas en este else, por lo que se conviente en progresivamente mas lento
            transcripts = set(row[0].split(','))
            if transcript not in transcripts:
                transcripts.add(transcript)
                updated = ",".join(transcripts)
                cursor.execute(
                    """
                    UPDATE kmers
                    SET transcripts = ?
                    WHERE kmer = ?
                    """,
                    (updated, kmer)
                )

def kmerize(seq, transcript, kmers, ksize, offset):
    seq = seq.upper()
    length = len(seq)

    for i in range(length - site_length + 1):
        window = seq[i:i + site_length]
        kmer = window[offset:offset + ksize]

        if kmer not in kmers:
            kmers[kmer] = set()

        kmers[kmer].add(transcript)

def write_output(kmers, fasta, ksize):
    outfile = fasta.parent / f"{fasta.name}.{ksize}nt-kmers.tab"
    with open(outfile, "w") as out:
        for kmer, names in kmers.items():
            out.write(kmer + '\t' + ','.join(names) + '\n')
    
    return outfile

def load_annotation(cursor, filepath):
    """Load annotation TSV file into SQLite."""
    print("Loading annotation...")

    with open(filepath) as f:
        data = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            transcript, description = line.split('\t')
            data.append((transcript, description))

    cursor.executemany(
        "INSERT INTO annotation (transcript, description) VALUES (?, ?)",
        data
    )

    print(f"{len(data)} annotations loaded")

def process_fasta(fasta, ksize):
    offset = site_length - ksize - 1
    kmers = {}

    transcript = None
    seq = ""
    count = 0

    with open(fasta) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if transcript is not None:
                    kmerize(seq, transcript, kmers, ksize, offset)
                    count += 1
                    print("Processed transcripts for kmer database: " + str(count), end="\r", file=sys.stderr)

                transcript = line[1:].split()[0]
                seq = ""
            else:
                seq += line

        # process last sequence
        if transcript is not None:
            kmerize(seq, transcript, kmers, ksize, offset)
            count += 1

    print("Processed " + str(count), file=sys.stderr)
    return kmers

def add_species(fasta, ksize):
    kmers = process_fasta(fasta, ksize)
    outfile = write_output(kmers, fasta, ksize)
    return outfile

def create_kmer_db(db_file, fasta, ksize, annot_file):

    conn, cursor = connect_sqlite(db_file)
    create_tables(cursor)
    process_fasta_load(fasta, ksize, conn)
    load_annotation(cursor, annot_file)

    conn.commit()
    conn.close()


