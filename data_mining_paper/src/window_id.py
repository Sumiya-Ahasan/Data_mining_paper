"""Shared homology rules (decoy v2 + label v2), so decoys and labels use the SAME definition.

Presence rules for a reference gene vs a genome (BLAST HSPs):
  strict   : HSPs with pident >= 90 on one subject cover >= 80% of the gene (union of query intervals)
  fragment : not strict, but some HSP has a 150-column window with identity >= 0.90
             (a read-length piece of the gene is really there) -> excluded from scoring
  null     : everything else (distant homologs included: they are the hard, realistic nulls)
Decoys v2 are built so that every 150-bp window vs their parent has identity < 0.90, i.e. a decoy is
exactly "null" w.r.t. its parent by the rule above.
"""
import numpy as np

WIN = 150
FRAG_ID = 0.90
STRICT_ID, STRICT_COV = 90.0, 0.80


def max_window_identity(qseq, sseq, w=WIN):
    """Max identity over all w-column windows of an alignment (gap = mismatch). None if shorter than w."""
    if len(qseq) < w:
        return None
    q = np.frombuffer(qseq.upper().encode(), dtype=np.uint8)
    s = np.frombuffer(sseq.upper().encode(), dtype=np.uint8)
    m = ((q == s) & (q != ord("-"))).astype(np.int32)
    c = np.concatenate([[0], np.cumsum(m)])
    return float((c[w:] - c[:-w]).max() / w)


def union_len(intervals):
    """Total length of the union of [start, end] (1-based, inclusive) intervals."""
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None: tot += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None: tot += cur_e - cur_s + 1
    return tot


def classify(hsps):
    """hsps: list of dicts (sseqid, pident, qlen, qstart, qend, maxwin) for ONE query.
    Returns (status, best_cov90, best_maxwin) with status in {'strict', 'fragment', None}."""
    by_s, best_win = {}, 0.0
    for h in hsps:
        if h["pident"] >= STRICT_ID:
            by_s.setdefault(h["sseqid"], []).append((h["qstart"], h["qend"]))
        if h["maxwin"] is not None: best_win = max(best_win, h["maxwin"])
    qlen = hsps[0]["qlen"] if hsps else 1
    cov = max((union_len(v) / qlen for v in by_s.values()), default=0.0)
    if cov >= STRICT_COV: return "strict", round(cov, 3), round(best_win, 4)
    if best_win >= FRAG_ID: return "fragment", round(cov, 3), round(best_win, 4)
    return None, round(cov, 3), round(best_win, 4)
