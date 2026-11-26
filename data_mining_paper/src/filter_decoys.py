#!/usr/bin/env python3
"""Drop decoys that look like PRESENCE of any real CARD gene under the label v2 rule (window_id.py):
any HSP with a 150-column window at >= 90% identity, or strict coverage, vs ANY CARD homolog member.

Input: BLAST outfmt '6 qseqid sseqid pident length qlen qstart qend qseq sseq' (query = decoys).
Usage: python filter_decoys.py --blast decoys_vs_card.tsv --in decoys_novel.fasta --out decoys_v2.fasta
"""
import argparse, os, sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from window_id import max_window_identity, classify  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", required=True); ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    hs = defaultdict(list)
    for line in open(a.blast):
        q, s, pid, L, qlen, qs, qe, qseq, sseq = line.rstrip("\n").split("\t")
        hs[q].append(dict(sseqid=s, pident=float(pid), qlen=int(qlen), qstart=int(qs), qend=int(qe),
                          maxwin=max_window_identity(qseq, sseq)))
    bad = {q for q, h in hs.items() if classify(h)[0] is not None}
    kept, drop, name, keep = Counter(), Counter(), None, False
    with open(a.out, "w") as out:
        for line in open(a.inp):
            if line.startswith(">"):
                name = line[1:].split()[0]; keep = name not in bad
                t = name.split("|")[1][0]
                (kept if keep else drop)[t] += 1
            if keep: out.write(line)
    print("decoys kept by type:", dict(kept), "| dropped (look like a real CARD gene):", dict(drop))


if __name__ == "__main__":
    main()
