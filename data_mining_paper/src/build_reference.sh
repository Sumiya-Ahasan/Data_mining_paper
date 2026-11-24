#!/usr/bin/env bash
# Build card_90 targets + decoys and index for KMA and Bowtie2.
# Runs on Kaggle (Ubuntu, Internet ON) or any Linux box. ~10-20 min on 4 cores.
# Usage: bash build_reference.sh <outdir>   (repo checked out so src/ is next to this script)
set -euo pipefail
SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO=sudo   # Kaggle runs as root (no sudo)
OUT=${1:-ref}; SRC=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT" && cd "$OUT"

# 0. tools (apt has cd-hit, bowtie2, seqtk; KMA built from source)
if ! command -v cd-hit-est >/dev/null; then
  $SUDO apt-get -qq update && $SUDO apt-get -qq install -y cd-hit bowtie2 seqtk >/dev/null
fi
if ! command -v kma >/dev/null; then
  git clone -q https://bitbucket.org/genomicepidemiology/kma.git && make -s -C kma && export PATH="$PWD/kma:$PATH"
fi

# 1. CARD data (record the version you get in docs/decision_log.md)
if [ ! -f card/nucleotide_fasta_protein_homolog_model.fasta ]; then
  mkdir -p card && curl -sL https://card.mcmaster.ca/latest/data -o card/card-data.tar.bz2
  tar -xjf card/card-data.tar.bz2 -C card
fi
grep -m1 -o '"_version": *"[^"]*"' card/card.json > card_version.txt || true
echo "CARD: $(cat card_version.txt)"; grep -c '>' card/nucleotide_fasta_protein_homolog_model.fasta

# 2. card_90: cluster at 90% nucleotide identity
cd-hit-est -i card/nucleotide_fasta_protein_homolog_model.fasta -o card90.fasta \
  -c 0.90 -n 8 -d 0 -M 8000 -T 4 > cdhit_card90.log
echo "card_90 clusters: $(grep -c '>' card90.fasta)"   # expect ~1318 (doc value; depends on CARD version)

# 3. decoys (d=3 per target, 75-85% identity)
python3 "$SRC/make_decoys.py" --in card90.fasta --out decoys_raw.fasta --per-target 3 --seed 2026

# 4. drop decoys >=90% identical to ANY target, then dedupe decoys among themselves at 90%
cd-hit-est-2d -i card90.fasta -i2 decoys_raw.fasta -o decoys_novel.fasta \
  -c 0.90 -n 8 -d 0 -M 8000 -T 4 > cdhit_2d.log
cd-hit-est -i decoys_novel.fasta -o decoys.fasta -c 0.90 -n 8 -d 0 -M 8000 -T 4 > cdhit_decoys.log
echo "decoys raw: $(grep -c '>' decoys_raw.fasta)  after 2d: $(grep -c '>' decoys_novel.fasta)  final: $(grep -c '>' decoys.fasta)"

# 5. combined index
cat card90.fasta decoys.fasta > card90_decoys.fasta
mkdir -p idx
kma index -i card90_decoys.fasta -o idx/card90d
bowtie2-build --threads 4 -q card90_decoys.fasta idx/card90d
# targets-only index too (for the decoy-competition side-effect check in the plan)
kma index -i card90.fasta -o idx/card90
bowtie2-build --threads 4 -q card90.fasta idx/card90
echo "done -> $OUT/idx"
