#!/usr/bin/env bash
# Descarga las estructuras experimentales de referencia desde el Protein Data Bank.
set -euo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$DIR"

echo "Descargando Chignolin (CLN025, PDB 5AWL)..."
curl -fsSL "https://files.rcsb.org/download/5AWL.pdb" -o "$DIR/5awl.pdb"
echo "  -> $DIR/5awl.pdb"

# Objetivo ambicioso (descomenta cuando llegues a la fase Trp-cage):
# echo "Descargando Trp-cage (PDB 1L2Y)..."
# curl -fsSL "https://files.rcsb.org/download/1L2Y.pdb" -o "$DIR/1l2y.pdb"
# echo "  -> $DIR/1l2y.pdb"

echo "Listo."
