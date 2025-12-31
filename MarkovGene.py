#!/usr/bin/env python3
import os
import math
import csv
import argparse
from collections import defaultdict, Counter

DNA = set("ACGT")


def list_fna_files(folder: str):
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Folder not found: {folder}")
    out = []
    for name in sorted(os.listdir(folder)):
        if name.startswith("."):
            continue
        if not name.lower().endswith(".fna"):
            continue
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            out.append(path)
    return out


def read_fna(path: str) -> str:
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip().upper()
            if not line or line.startswith(">"):
                continue
            chunks.append("".join(c for c in line if c in DNA))
    return "".join(chunks)


def gene_name(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


class MarkovK:
    def __init__(self, k: int = 3, alpha: float = 1.0):
        self.k = k
        self.alpha = alpha
        self.counts = defaultdict(Counter)  # context -> next base counts

    def fit(self, sequences):
        pad = "^" * self.k
        for seq in sequences:
            if not seq:
                continue
            s = pad + seq
            for i in range(self.k, len(s)):
                ctx = s[i - self.k : i]
                nxt = s[i]
                if nxt in DNA:
                    self.counts[ctx][nxt] += 1

    def log_likelihood(self, seq: str) -> float:
        if not seq:
            return float("-inf")
        pad = "^" * self.k
        s = pad + seq
        ll = 0.0
        for i in range(self.k, len(s)):
            ctx = s[i - self.k : i]
            nxt = s[i]
            row = self.counts.get(ctx, Counter())
            total = sum(row.values())
            denom = total + self.alpha * 4
            p = (row.get(nxt, 0) + self.alpha) / denom
            ll += math.log(p)
        return ll


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_pos", default="train/pos")
    ap.add_argument("--train_neg", default="train/neg")
    ap.add_argument("--test_pos", default="test/pos")
    ap.add_argument("--test_neg", default="test/neg")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--out", default="scores.csv")
    args = ap.parse_args()

    train_pos_files = list_fna_files(args.train_pos)
    train_neg_files = list_fna_files(args.train_neg)
    if not train_pos_files or not train_neg_files:
        raise RuntimeError("Both train/pos and train/neg must contain at least one .fna file.")

    pos_seqs = [read_fna(p) for p in train_pos_files]
    neg_seqs = [read_fna(p) for p in train_neg_files]

    pos_model = MarkovK(k=args.k, alpha=args.alpha)
    neg_model = MarkovK(k=args.k, alpha=args.alpha)
    pos_model.fit(pos_seqs)
    neg_model.fit(neg_seqs)

    test_pos_files = list_fna_files(args.test_pos)
    test_neg_files = list_fna_files(args.test_neg)

    rows = []
    idx = 1
    for label, files in (("pos", test_pos_files), ("neg", test_neg_files)):
        for path in files:
            seq = read_fna(path)
            score = pos_model.log_likelihood(seq) - neg_model.log_likelihood(seq)
            rows.append((idx, label, gene_name(path), score))
            idx += 1

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["number", "label", "gene_name", "score"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
