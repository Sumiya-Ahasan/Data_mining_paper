"""Validate subsampled paired FASTQ.gz inside a zip: counts, pairing, format, lengths, spot spread, GC, N, quality."""
import gzip, io, sys, zipfile, re, json, statistics

zp, manifest, out = sys.argv[1], sys.argv[2], sys.argv[3]
nreads = {}
with open(manifest, encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t"); nreads[f[0]] = int(f[4])

z = zipfile.ZipFile(zp)
runs = sorted({re.sub(r"_[12]\.fastq\.gz$", "", n.split("/")[-1]) for n in z.namelist() if n.endswith(".fastq.gz")})
report = {}
for run in runs:
    f1 = gzip.open(io.BufferedReader(z.open(f"bg/{run}_1.fastq.gz"), 1 << 20), "rt")
    f2 = gzip.open(io.BufferedReader(z.open(f"bg/{run}_2.fastq.gz"), 1 << 20), "rt")
    n = mism = badfmt = 0; lens1 = {}; lens2 = {}; gc = at = nn = 0; qsum = qn = 0
    spots = []; dup = 0; seen_prev = None; nonmono = 0
    while True:
        h1 = f1.readline()
        h2 = f2.readline()
        if not h1 or not h2:
            extra = bool(h1) or bool(h2); break
        s1, p1, q1 = f1.readline().rstrip(), f1.readline(), f1.readline().rstrip()
        s2, p2, q2 = f2.readline().rstrip(), f2.readline(), f2.readline().rstrip()
        n += 1
        id1 = h1.split()[0].rstrip("\n"); id2 = h2.split()[0].rstrip("\n")
        id1 = re.sub(r"/1$", "", id1); id2 = re.sub(r"/2$", "", id2)
        if id1 != id2: mism += 1
        if not (h1.startswith("@") and p1.startswith("+") and len(s1) == len(q1) and h2.startswith("@") and p2.startswith("+") and len(s2) == len(q2)):
            badfmt += 1
        lens1[len(s1)] = lens1.get(len(s1), 0) + 1; lens2[len(s2)] = lens2.get(len(s2), 0) + 1
        m = re.search(r"\.(\d+)$", id1)
        if m:
            sp = int(m.group(1))
            if seen_prev is not None and sp <= seen_prev: nonmono += 1
            seen_prev = sp
            if n % 1000 == 1: spots.append(sp)
        if n <= 200000:
            for s in (s1, s2):
                gc += s.count("G") + s.count("C"); at += s.count("A") + s.count("T"); nn += s.count("N")
            for q in (q1, q2):
                qsum += sum(ord(c) - 33 for c in q); qn += len(q)
    top = lambda d: sorted(d.items(), key=lambda x: -x[1])[:3]
    tot = nreads.get(run)
    r = dict(pairs=n, unequal_file_length=extra, name_mismatch=mism, bad_format=badfmt,
             len_r1_top=top(lens1), len_r2_top=top(lens2),
             gc_pct=round(100 * gc / max(1, gc + at), 2), n_pct=round(100 * nn / max(1, gc + at + nn), 3),
             mean_q=round(qsum / max(1, qn), 1),
             spot_min=min(spots) if spots else None, spot_max=max(spots) if spots else None,
             total_pairs_in_run=tot, nonmonotonic_ids=nonmono,
             spot_span_frac=round((max(spots) - min(spots)) / tot, 3) if spots and tot else None)
    report[run] = r
    print(run, json.dumps(r), flush=True)
json.dump(report, open(out, "w"), indent=1)
