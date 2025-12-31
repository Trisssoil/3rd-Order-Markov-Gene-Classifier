#!/usr/bin/env python3
import random
from pathlib import Path
from typing import List, Tuple

BASES = set("ACGT")

# ---------- FASTA I/O ----------
def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records = []
    name = None
    chunks = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if name is not None:
            records.append((name, "".join(chunks)))
    return records

def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def clean_seq(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in BASES)

# ---------- Shuffling ----------
def shuffle_sequence(seq: str, rng: random.Random) -> str:
    """
    Mononucleotide-preserving shuffle.
    Keeps base composition and length.
    """
    seq = clean_seq(seq)
    if len(seq) <= 1:
        return seq
    arr = list(seq)
    rng.shuffle(arr)
    return "".join(arr)

# ---------- Main ----------
def main():
    base = Path(__file__).resolve().parent
    coding_dir = base / "folder_coding"
    noncoding_dir = base / "folder_nonCoding"

    if not coding_dir.exists():
        raise FileNotFoundError(f"Missing folder: {coding_dir}")

    noncoding_dir.mkdir(exist_ok=True)

    fasta_exts = {".fa", ".fasta", ".fna", ".ffn", ".fas"}
    coding_files = sorted(
        p for p in coding_dir.iterdir()
        if p.is_file() and p.suffix.lower() in fasta_exts
    )

    if not coding_files:
        raise RuntimeError(f"No FASTA files found in {coding_dir}")

    rng = random.Random(42)  # fixed seed for reproducibility

    for infile in coding_files:
        records = read_fasta(infile)
        shuffled_records = []

        for name, seq in records:
            shuffled_seq = shuffle_sequence(seq, rng)
            shuffled_records.append((name, shuffled_seq))

        # 🔑 change: prefix filename with "shuffled_"
        outfile = noncoding_dir / f"shuffled_{infile.name}"
        write_fasta(outfile, shuffled_records)

        print(f"Shuffled: {infile.name} -> folder_nonCoding/{outfile.name}")

    print("\nDone. Shuffled noncoding files created.")

if __name__ == "__main__":
    main()
