#!/usr/bin/env python3
"""
Create a SQLite database and load TSV data into tables.
"""

import sqlite3
import sys


# -------------------------------
# Functions
# -------------------------------

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


def load_kmers(cursor, filepath, batch_size=50000):
    """Load kmers TSV file into SQLite using streaming batches."""

    print("Loading kmers (streaming)...")

    batch = []
    count = 0

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            kmer, transcripts = line.split('\t')
            batch.append((kmer, transcripts))

            if len(batch) >= batch_size:
                cursor.executemany(
                    "INSERT INTO kmers (kmer, transcripts) VALUES (?, ?)",
                    batch
                )
                count += len(batch)
                print(f"\r{count} kmers loaded...", end="", flush=True)
                batch.clear()

        # flush last batch
        if batch:
            cursor.executemany(
                "INSERT INTO kmers (kmer, transcripts) VALUES (?, ?)",
                batch
            )
            count += len(batch)

    print(f"Done: {count} kmers loaded")


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

def create_sqlite_db(db_file, kmers_file, annotation_file):
    conn, cursor = connect_sqlite(db_file)

    create_tables(cursor)

    load_kmers(cursor, kmers_file)
    load_annotation(cursor, annotation_file)

    conn.commit()
    conn.close()

    print("DB created")

