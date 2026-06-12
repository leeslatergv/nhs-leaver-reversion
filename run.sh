#!/usr/bin/env bash
# Run the whole pipeline, from the project root: ./run.sh
# The scripts are numbered in dependency order, so this just runs them in turn.
set -euo pipefail
cd "$(dirname "$0")"

for f in R/0[1-9]_*.R; do
  echo "==> $f"
  Rscript "$f"
done
