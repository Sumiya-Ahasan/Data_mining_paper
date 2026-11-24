#!/usr/bin/env python3
"""Generate homology-preserving decoy references for card_90 cluster representatives.

Each decoy is a substitution-only mutant of one target representative, with
nucleotide identity drawn uniformly from [--min-id, --max-id] (default 0.75-0.85).
Substitutions are biased to 3rd codon positions and never create an in-frame stop
codon, so decoys look like plausible divergent homologs, not random sequence.

Decoys are null by construction. Run check_decoys.sh afterwards to drop any decoy
that is >= 90% identical to any target or to another decoy (cd-hit-est-2d / cd-hit-est).

Usage:
  python make_decoys.py --in card90_rep.fasta --out decoys_raw.fasta --per-target 3 --seed 2026
"""
import argparse
import random

STOPS = {"TAA", "TAG", "TGA"}
BASES = "ACGT"


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def creates_stop(seq, pos, new_base):
    """True if placing new_base at pos creates a stop codon in frame 0."""
    start = pos - pos % 3
    codon = seq[start:start + 3]
    if len(codon) < 3:
        return False
    codon[pos - start] = new_base
    return "".join(codon) in STOPS


def mutate(seq, identity, rng, third_pos_weight=0.7):
    s = list(seq)
    n = len(s)
    n_sub = round(n * (1.0 - identity))
    # positions eligible for substitution: only A/C/G/T
    pos3 = [i for i in range(n) if i % 3 == 2 and s[i] in BASES]
    pos12 = [i for i in range(n) if i % 3 != 2 and s[i] in BASES]
    rng.shuffle(pos3)
    rng.shuffle(pos12)
    done, guard = set(), 0
    while len(done) < n_sub and guard < 50 * n:
        guard += 1
        pool = pos3 if (rng.random() < third_pos_weight and pos3) else pos12
        if not pool:
            pool = pos3 or pos12
        if not pool:
            break
        i = pool.pop()
        alts = [b for b in BASES if b != s[i] and not creates_stop(s, i, b)]
        if not alts:
            continue
        s[i] = rng.choice(alts)
        done.add(i)
    realised = 1.0 - len(done) / n
    return "".join(s), realised


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-target", type=int, default=3)
    ap.add_argument("--min-id", type=float, default=0.75)
    ap.add_argument("--max-id", type=float, default=0.85)
    ap.add_argument("--min-len", type=int, default=150, help="skip targets shorter than this")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    n_t = n_d = 0
    with open(a.out, "w") as out, open(a.out + ".tsv", "w") as log:
        log.write("decoy_id\ttarget_id\ttarget_len\ttarget_identity\n")
        for name, seq in read_fasta(a.inp):
            n_t += 1
            if len(seq) < a.min_len:
                continue
            for k in range(a.per_target):
                ident = rng.uniform(a.min_id, a.max_id)
                dseq, real = mutate(seq, ident, rng)
                did = f"DECOY|{k}|{name}"
                out.write(f">{did}\n")
                for j in range(0, len(dseq), 80):
                    out.write(dseq[j:j + 80] + "\n")
                log.write(f"{did}\t{name}\t{len(seq)}\t{real:.4f}\n")
                n_d += 1
    print(f"targets read: {n_t}; decoys written: {n_d} -> {a.out}")


if __name__ == "__main__":
    main()
