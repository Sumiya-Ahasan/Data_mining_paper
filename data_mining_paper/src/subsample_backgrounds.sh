#!/usr/bin/env bash
# Download each background run from ENA (resumable, retried), keep exactly N read pairs, delete raw.
# Raw files go to $TMPDIR (default /tmp, large on Kaggle); outputs to <outdir>. Same seed on R1/R2 keeps pairs.
#
# Usage: bash subsample_backgrounds.sh <manifest.tsv> <outdir> [N] [HEAD_PAIRS]
#   manifest: ENA filereport TSV, col1 = run_accession, col5 = read_count (pairs)
#   HEAD_PAIRS (optional): use only the first HEAD_PAIRS pairs of each run (huge PRJNA1095031 runs,
#     e.g. 50000000). Those are streamed (no full download); a stream is retried up to 3 times.
set -uo pipefail
SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO=sudo   # Kaggle runs as root (no sudo)
MAN=$1; OUT=$2; N=${3:-5000000}; HP=${4:-0}
RAW=${TMPDIR:-/tmp}/raw_bg; mkdir -p "$OUT" "$RAW"
command -v seqtk >/dev/null || { $SUDO apt-get -qq update && $SUDO apt-get -qq install -y seqtk >/dev/null; }

lines() { zcat "$1" 2>/dev/null | wc -l; }

tail -n +2 "$MAN" | cut -f1,5 | while read -r RUN NREADS; do
  O1="$OUT/${RUN}_1.fastq.gz"; O2="$OUT/${RUN}_2.fastq.gz"
  if [ -s "$O1" ] && [ "$(lines "$O1")" -eq $((4*N)) ] && [ "$(lines "$O2")" -eq $((4*N)) ]; then
    echo "skip $RUN (complete)"; continue; fi
  URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2)
  U1=$(echo "$URLS" | cut -d';' -f1); U2=$(echo "$URLS" | cut -d';' -f2)

  for attempt in 1 2 3; do
    if [ "$HP" -gt 0 ] && [ "$HP" -lt "$NREADS" ]; then
      POOL=$HP; FRAC=$(python3 -c "print(min(1.0, 1.03*$N/$POOL))")
      echo "$RUN stream first $POOL pairs, frac=$FRAC (attempt $attempt)"
      for r in 1 2; do U=$([ $r = 1 ] && echo "$U1" || echo "$U2")
        ( curl -sS --retry 5 "https://$U" | gzip -dc 2>/dev/null | head -n $((4*POOL)) \
          | seqtk sample -s 11 - "$FRAC" | head -n $((4*N)) | gzip -1 > "$OUT/${RUN}_$r.fastq.gz" ) &
      done; wait
    else
      FRAC=$(python3 -c "print(min(1.0, 1.03*$N/$NREADS))")
      echo "$RUN download full files, frac=$FRAC (attempt $attempt)"
      for r in 1 2; do U=$([ $r = 1 ] && echo "$U1" || echo "$U2")
        ( wget -q -c --tries=20 --timeout=60 -O "$RAW/${RUN}_$r.fastq.gz" "https://$U" ) &
      done; wait
      # integrity check of downloads before sampling
      for r in 1 2; do gzip -t "$RAW/${RUN}_$r.fastq.gz" || echo "  corrupt download R$r"; done
      for r in 1 2; do
        ( seqtk sample -s 11 "$RAW/${RUN}_$r.fastq.gz" "$FRAC" | head -n $((4*N)) | gzip -1 > "$OUT/${RUN}_$r.fastq.gz" ) &
      done; wait
    fi
    L1=$(lines "$O1"); L2=$(lines "$O2")
    echo "  R1 lines: $L1  R2 lines: $L2 (expect $((4*N)) each)"
    if [ "$L1" -eq $((4*N)) ] && [ "$L2" -eq $((4*N)) ]; then
      bad=$(paste <(zcat "$O1" | awk 'NR%4==1{sub(/\/1$/,"",$1);print $1}') \
                  <(zcat "$O2" | awk 'NR%4==1{sub(/\/2$/,"",$1);print $1}') | awk '$1!=$2{b++} END{print b+0}')
      echo "  name mismatches (all pairs): $bad (expect 0)"
      [ "$bad" -eq 0 ] && { rm -f "$RAW/${RUN}"_*.fastq.gz; echo "  OK $RUN"; break; }
    fi
    echo "  FAILED attempt $attempt for $RUN"; rm -f "$O1" "$O2"
    [ "$attempt" -eq 3 ] && echo "  GIVING UP on $RUN — rerun later"
  done
done
echo "ALL DONE"; ls -la "$OUT"
