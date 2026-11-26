#!/usr/bin/env bash
# Reference v2 = SAME CARD 4.0.2 card_90 targets as v1 + decoy v2 (types U/N/M, see make_decoys_v2.py).
# Reuses card/, card90.fasta and card90.fasta.clstr from the v1 reference so targets do not change.
# Usage: bash build_reference_v2.sh <v1_ref_dir> <outdir>      (needs cd-hit, ncbi-blast+, bowtie2, $KMA or kma)
set -euo pipefail
V1=$1; OUT=$2; SRC=$(cd "$(dirname "$0")" && pwd); KMA=${KMA:-kma}
mkdir -p "$OUT" && cd "$OUT"
cp -r "$V1/card" . && cp "$V1/card90.fasta" "$V1/card90.fasta.clstr" .
[ -f "$V1/card_version.txt" ] && cp "$V1/card_version.txt" . || true
echo "targets: $(grep -c '>' card90.fasta)"

# 1. decoys v2: per target 1 uniform (v1-like), 2 near, 2 mosaic
python3 "$SRC/make_decoys_v2.py" --in card90.fasta --out decoys_raw.fasta --n-U 1 --n-N 2 --n-M 2 --seed 2026

# 2. drop decoys >= 90% to any target representative (cd-hit-est-2d, as v1)
cd-hit-est-2d -i card90.fasta -i2 decoys_raw.fasta -o decoys_novel.fasta -c 0.90 -n 8 -d 0 -M 8000 -T 4 > cdhit_2d.log

# 3. drop decoys that would count as present for ANY of the 6059 CARD members (label v2 rule, 150-bp windows)
makeblastdb -in card/nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out card_members > /dev/null
blastn -task dc-megablast -query decoys_novel.fasta -db card_members -evalue 1e-10 -max_target_seqs 50 -num_threads 4 \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend qseq sseq' > decoys_vs_card.tsv
python3 "$SRC/filter_decoys.py" --blast decoys_vs_card.tsv --in decoys_novel.fasta --out decoys_filtered.fasta
rm -f decoys_vs_card.tsv card_members.*

# 4. dedupe decoys among themselves at 90%
cd-hit-est -i decoys_filtered.fasta -o decoys.fasta -c 0.90 -n 8 -d 0 -M 8000 -T 4 > cdhit_decoys.log
echo "decoys raw: $(grep -c '>' decoys_raw.fasta)  after 2d: $(grep -c '>' decoys_novel.fasta)  after CARD filter: $(grep -c '>' decoys_filtered.fasta)  final: $(grep -c '>' decoys.fasta)"
for t in U N M; do echo "  type $t: $(grep -c "^>DECOY|$t" decoys.fasta)"; done

# 5. combined index (same file names as v1, so pilot.py/features.py work unchanged with --ref <outdir>)
cat card90.fasta decoys.fasta > card90_decoys.fasta
mkdir -p idx
$KMA index -i card90_decoys.fasta -o idx/card90d
bowtie2-build --threads 4 -q card90_decoys.fasta idx/card90d
echo "REF V2 DONE -> $OUT"
