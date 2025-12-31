#!/usr/bin/env python3

import os
import random

INPUT_DIR = "Genes"
OUTPUT_DIR = "Controls"


def read_sequence(filepath):
    """Read a gene sequence, ignoring FASTA headers if present."""
    seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq)


def shuffle_sequence(seq):
    """Shuffle a DNA sequence while preserving nucleotide composition."""
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)


def write_sequence(filepath, seq):
    """Write sequence in plain text (single line)."""
    with open(filepath, "w") as f:
        f.write(seq + "\n")


def main():
    if not os.path.exists(INPUT_DIR):
        raise FileNotFoundError(f"Input folder '{INPUT_DIR}' not found.")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    gene_files = sorted(
        f for f in os.listdir(INPUT_DIR)
        if os.path.isfile(os.path.join(INPUT_DIR, f))
    )

    if not gene_files:
        print("No gene files found to shuffle.")
        return

    for filename in gene_files:
        input_path = os.path.join(INPUT_DIR, filename)

        seq = read_sequence(input_path)
        if not seq:
            print(f"Skipping empty file: {filename}")
            continue

        shuffled_seq = shuffle_sequence(seq)

        output_filename = f"shuffled_{filename}"
        output_path = os.path.join(OUTPUT_DIR, output_filename)

        write_sequence(output_path, shuffled_seq)

        print(f"Shuffled: {filename} → {output_filename}")

    print("All test genes shuffled successfully.")


if __name__ == "__main__":
    main()
