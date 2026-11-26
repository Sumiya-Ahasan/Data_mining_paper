#!/usr/bin/env python3
"""Decoy v2: decoys that mimic how real absent-ARG homologs look (G2 extension, 2026-11-26).

The G2 extension showed uniform 75-85% decoys are too easy: absent ARGs get reads from real homologs that
have conserved blocks (reads map at 85-90% identity locally). v2 makes three decoy types per target:
  U  uniform   : v1-style, identity U(0.75, 0.85), random positions       (kept for comparison)
  N  near      : identity U(0.86, 0.89), mutations spread evenly         (cluster neighbours just under 90%)
  M  mosaic    : conserved blocks (identity 0.87-0.89, 150-400 bp) alternating with variable blocks
                 (identity 0.55-0.70, 100-400 bp)                         (distant homologs with conserved domains)
All types: substitutions only, 3rd codon positions preferred, never an in-frame stop, and every 150-bp
window vs the parent has identity < 0.90 (window_id.FRAG_ID) -> a decoy is "null" by the label v2 rule.
Ids: DECOY|<type><k>|<target id>.  A .tsv log gives type, global identity and max window identity.

Usage: python make_decoys_v2.py --in card90.fasta --out decoys_raw.fasta --n-U 1 --n-N 2 --n-M 2 --seed 2026
"""
import argparse, os, random, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from window_id import WIN, FRAG_ID  # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
BASES = "ACGT"
MIN_MM = round(WIN * (1 - FRAG_ID)) + 1       # >= 16 mismatches per 150 bp  ->  window identity <= 0.893


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq).upper()
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name is not None: yield name, "".join(seq).upper()


def try_sub(s, i, rng):
    """Substitute s[i] (list) with a base that differs and creates no in-frame stop. True if done."""
    if s[i] not in BASES: return False
    st = i - i % 3
    alts = []
    for b in BASES:
        if b == s[i]: continue
        cod = s[st:st + 3]
        if len(cod) == 3:
            cod = cod.copy(); cod[i - st] = b
            if "".join(cod) in STOPS: continue
        alts.append(b)
    if not alts: return False
    s[i] = rng.choice(alts); return True


def mutate_region(s, mut, start, end, rate, rng, spread, w3=0.7):
    """Mutate about rate*(end-start) positions of [start, end). spread=True: one per slot of length 1/rate."""
    n = round((end - start) * rate)
    if n <= 0: return
    if spread:
        step = (end - start) / n
        for k in range(n):
            lo, hi = start + int(k * step), max(start + int(k * step) + 1, start + int((k + 1) * step))
            cand = [i for i in range(lo, min(hi, end)) if not mut[i]]
            c3 = [i for i in cand if i % 3 == 2]
            rng.shuffle(cand)
            order = (c3 if (c3 and rng.random() < w3) else []) + cand
            for i in order:
                if try_sub(s, i, rng): mut[i] = True; break
    else:
        cand = [i for i in range(start, end) if not mut[i]]
        p3 = [i for i in cand if i % 3 == 2]; p12 = [i for i in cand if i % 3 != 2]
        rng.shuffle(p3); rng.shuffle(p12)
        done = 0
        while done < n and (p3 or p12):
            pool = p3 if (p3 and (rng.random() < w3 or not p12)) else p12
            i = pool.pop()
            if try_sub(s, i, rng): mut[i] = True; done += 1


def cap_windows(s, mut, rng):
    """Add mutations until every WIN-bp window has >= MIN_MM mismatches vs the parent."""
    L = len(s)
    if L < WIN: return
    for _ in range(10 * L):
        c = [0]
        for x in mut: c.append(c[-1] + x)
        bad = next((j for j in range(L - WIN + 1) if c[j + WIN] - c[j] < MIN_MM), None)
        if bad is None: return
        cand = [i for i in range(bad, bad + WIN) if not mut[i]]
        rng.shuffle(cand)
        for i in sorted(cand, key=lambda i: i % 3 != 2)[:4]:   # 3rd positions first
            if try_sub(s, i, rng): mut[i] = True; break
        else:
            return


def max_win_id(mut):
    L = len(mut)
    if L < WIN: return None
    c = [0]
    for x in mut: c.append(c[-1] + x)
    return 1 - min(c[j + WIN] - c[j] for j in range(L - WIN + 1)) / WIN


def make(seq, kind, rng):
    s, mut = list(seq), [False] * len(seq)
    L = len(seq)
    if kind == "U":
        mutate_region(s, mut, 0, L, 1 - rng.uniform(0.75, 0.85), rng, spread=False)
    elif kind == "N":
        mutate_region(s, mut, 0, L, 1 - rng.uniform(0.86, 0.89), rng, spread=True)
    elif kind == "M":
        pos, conserved = 0, rng.random() < 0.5
        while pos < L:
            blen = rng.randint(150, 400) if conserved else rng.randint(100, 400)
            end = min(L, pos + blen)
            if conserved: mutate_region(s, mut, pos, end, 1 - rng.uniform(0.87, 0.89), rng, spread=True)
            else: mutate_region(s, mut, pos, end, 1 - rng.uniform(0.55, 0.70), rng, spread=False)
            pos, conserved = end, not conserved
    cap_windows(s, mut, rng)
    return "".join(s), 1 - sum(mut) / L, max_win_id(mut)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-U", type=int, default=1)
    ap.add_argument("--n-N", type=int, default=2)
    ap.add_argument("--n-M", type=int, default=2)
    ap.add_argument("--min-len", type=int, default=150, help="skip targets shorter than this")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    n_t = n_d = 0
    with open(a.out, "w") as out, open(a.out + ".tsv", "w") as log:
        log.write("decoy_id\ttarget_id\ttype\ttarget_len\tidentity\tmax_window_identity\n")
        for name, seq in read_fasta(a.inp):
            n_t += 1
            if len(seq) < a.min_len: continue
            for kind, n in (("U", a.n_U), ("N", a.n_N), ("M", a.n_M)):
                for k in range(n):
                    d, ident, mw = make(seq, kind, rng)
                    did = f"DECOY|{kind}{k}|{name}"
                    out.write(f">{did}\n" + "".join(d[j:j + 80] + "\n" for j in range(0, len(d), 80)))
                    log.write(f"{did}\t{name}\t{kind}\t{len(seq)}\t{ident:.4f}\t{mw:.4f}\n")
                    n_d += 1
    print(f"targets read: {n_t}; decoys written: {n_d} -> {a.out}")


if __name__ == "__main__":
    main()
