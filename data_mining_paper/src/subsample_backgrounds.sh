#!/usr/bin/env bash
# Stream each background run from ENA and keep exactly N read pairs (default 5M),
# without storing the full FASTQ (Kaggle disk is limited). Same seed on R1/R2 keeps pairs.
#
# Usage: bash subsample_backgrounds.sh <manifest.tsv> <outdir> [N] [HEAD_PAIRS]
#   manifest: ENA filereport TSV, col1 = run_accession, col5 = read_count (pairs)
#   HEAD_PAIRS (optional): stream only the first HEAD_PAIRS pairs of each run, then sample N
#     at random from those. Use for the huge PRJNA1095031 runs (e.g. 50000000).
#     Known bias: reads at the start of a file come from the same flowcell tiles; report it.
# Examples:
#   bash src/subsample_backgrounds.sh data_manifest/bg_hiseq_prjeb13831.tsv bg 5000000
#   bash src/subsample_backgrounds.sh data_manifest/bg_test_prjna1095031.tsv bg 5000000 50000000
set -euo pipefail
SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO=sudo   # Kaggle runs as root (no sudo)
MAN=$1; OUT=$2; N=${3:-5000000}; HP=${4:-0}
mkdir -p "$OUT"
command -v seqtk >/dev/null || { $SUDO apt-get -qq update && $SUDO apt-get -qq install -y seqtk >/dev/null; }

tail -n +2 "$MAN" | cut -f1,5 | while read -r RUN NREADS; do
  [ -s "$OUT/${RUN}_2.fastq.gz" ] && { echo "skip $RUN"; continue; }
  POOL=$NREADS
  if [ "$HP" -gt 0 ] && [ "$HP" -lt "$NREADS" ]; then POOL=$HP; fi
  FRAC=$(python3 -c "print(min(1.0, 1.03*$N/$POOL))")   # 3% headroom, then truncate
  URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2)
  U1=$(echo "$URLS" | cut -d';' -f1); U2=$(echo "$URLS" | cut -d';' -f2)
  echo "$RUN reads=$NREADS pool=$POOL frac=$FRAC"
  # Identical seed and identical input order => identical read selection in R1 and R2.
  # '|| true': head closes pipes early by design (SIGPIPE); the checks below catch real failures.
  pull() {  # $1 = url, $2 = output file
    if [ "$POOL" -lt "$NREADS" ]; then
      curl -s "https://$1" | gzip -dc | head -n $((4*POOL)) | seqtk sample -s 11 - "$FRAC" | head -n $((4*N)) | gzip -1 > "$2"
    else
      curl -s "https://$1" | seqtk sample -s 11 - "$FRAC" | head -n $((4*N)) | gzip -1 > "$2"
    fi
  }
  { pull "$U1" "$OUT/${RUN}_1.fastq.gz" || true; } &
  { pull "$U2" "$OUT/${RUN}_2.fastq.gz" || true; }
  wait
  # sanity: both files must have 4*N lines and matching read names
  for r in 1 2; do echo "  R$r lines: $(zcat "$OUT/${RUN}_$r.fastq.gz" | wc -l) (expect $((4*N)))"; done
  paste <(zcat "$OUT/${RUN}_1.fastq.gz" | awk 'NR%4==1{print $1}' | head -1000 | sed 's#/1$##') \
        <(zcat "$OUT/${RUN}_2.fastq.gz" | awk 'NR%4==1{print $1}' | head -1000 | sed 's#/2$##') \
    | awk '$1!=$2{bad++} END{print "  first-1000 name mismatches: " bad+0 " (expect 0)"}'
done
