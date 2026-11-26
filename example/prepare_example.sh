#!/usr/bin/env bash
# Build the small OuroAssemble example: public complete genomes + simulated HiFi reads.
#   bash example/prepare_example.sh      (writes example/data/, ~50 MB)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p genomes data

# accession  depth   description
GENOMES=(
  "NC_000908.2 40"   # Mycoplasmoides genitalium G37 chromosome (580 kb)
  "NC_002528.1 30"   # Buchnera aphidicola APS chromosome (641 kb)
  "NC_002252.1 80"   # Buchnera aphidicola APS plasmid pLeu (7.3 kb)
  "NC_002253.1 80"   # Buchnera aphidicola APS plasmid pTrp (7.8 kb)
  "NC_001422.1 60"   # Escherichia phage phiX174 (5.4 kb, circular)
)

specs=()
for entry in "${GENOMES[@]}"; do
  read -r acc depth <<<"$entry"
  fa="genomes/${acc}.fa"
  if [ ! -s "$fa" ]; then
    curl -sSf --retry 3 -o "$fa" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
    sleep 0.5
  fi
  specs+=("${fa}:${depth}")
done

python3 simulate_hifi.py data/reads.fq.gz "${specs[@]}"
echo "Example reads: $(pwd)/data/reads.fq.gz"
