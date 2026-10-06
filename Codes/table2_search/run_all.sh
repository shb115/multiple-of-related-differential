#!/bin/bash
# Section 2.2: exhaustive check that Table 2 is complete for 4 <= w <= 8.
# Usage: bash run_all.sh [out_dir]      (default out_dir: the current directory)
# Builds lemma7_search in a temporary directory and writes results_mode0.txt and
# results_mode1.txt to out_dir; they match Results/table2_search/results_mode{0,1}.txt
# except for the timing lines.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$(mkdir -p "${1:-.}" && cd "${1:-.}" && pwd)"
BUILD="$(mktemp -d)"; trap 'rm -rf "$BUILD"' EXIT
gcc -O3 -march=native -o "$BUILD/lemma7_search" "$HERE/lemma7_search.c"
: > "$OUT/results_mode0.txt"; : > "$OUT/results_mode1.txt"
if [ -x /usr/bin/time ]; then
  run() { /usr/bin/time -f "wall %e s, maxrss %M KB" "$BUILD/lemma7_search" "$@"; }
else
  run() { "$BUILD/lemma7_search" "$@"; }
fi
# mode 0: exactly one zero per byte (the form of Table 2); raw = 1 adds an unnormalized cross-check for w <= 6
for spec in "4 0x13 1" "5 0x25 1" "6 0x43 1" "7 0x83 0" "8 0x11B 0"; do
  set -- $spec
  { run $1 $2 0 $3; } >> "$OUT/results_mode0.txt" 2>&1
done
# mode 1: at least one zero per byte (Definition 1)
for spec in "4 0x13" "5 0x25" "6 0x43" "7 0x83" "8 0x11B"; do
  set -- $spec
  { run $1 $2 1 0; } >> "$OUT/results_mode1.txt" 2>&1
done
