#!/usr/bin/env bash
# Run the T-BMT tests with one EPOCH2D binary.
#
#   ./run_tests.sh LABEL [EPOCH_BINARY] [DX_SCALE]
#
# Each test's input.deck is copied to <test>/run_LABEL/ and run there.
# DX_SCALE multiplies the cell size, and hence dt (CFL), for a dt scan.
# The binary's src/ is linked into each run directory so that the default
# table paths (qed_table_location etc.) resolve, as job_epoch does.
# Examples:
#   ./run_tests.sh hc                         # default binary (HC_PUSH)
#   ./run_tests.sh boris /path/to/boris/epoch2d
#   ./run_tests.sh hc_dx8 ../../../bin/epoch2d 8
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
label=${1:?usage: $0 LABEL [EPOCH_BINARY] [DX_SCALE]}
epoch=$(realpath "${2:-$here/../../../bin/epoch2d}")
scale=${3:-1}
nprocs=${MPIPROCS:-1}

for deck in "$here"/test*/input.deck; do
  tdir=$(dirname "$deck")
  rdir=$tdir/run_$label
  rm -rf "$rdir"
  mkdir -p "$rdir"
  sed "s/^  cell = \([0-9.]*\) \* micron/  cell = \1 * $scale * micron/" \
    "$deck" > "$rdir/input.deck"
  ln -s "$(dirname "$epoch")/../src" "$rdir/src"
  echo "== $(basename "$tdir") [$label, dx x $scale]"
  (cd "$rdir" && echo . | mpirun -n "$nprocs" "$epoch" > epoch.log 2>&1) \
    || { echo "   FAILED, see $rdir/epoch.log"; exit 1; }
  grep -h "Final runtime" "$rdir/epoch.log" || true
done
