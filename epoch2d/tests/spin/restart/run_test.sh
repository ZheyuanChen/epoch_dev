#!/usr/bin/env bash
# Run the spin restart test with one EPOCH2D binary.
#
#   ./run_test.sh LABEL [EPOCH_BINARY]
#
# Under run_LABEL/:
#   full/       input.deck on 2 ranks, 0-20 fs (dumps 0000-0010; 0005 is a
#               restart dump)
#   restart2/   restarted from full/0005.sdf on 2 ranks (same decomposition)
#   restart1/   restarted from full/0005.sdf on 1 rank (different
#               decomposition, so every particle is read onto another rank)
# Rank counts are kept at 1-2: EPOCH (upstream 4.21-devel too) writes
# corrupted field blocks into restart dumps at 3-4 ranks on this deck and
# then crashes or diverges on restart, with or without -DSPIN. See README.
# The binary's src/ is linked into each run directory so that the QED
# tables resolve, as job_epoch does.
# Then: python analyse.py LABEL
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
label=${1:?usage: $0 LABEL [EPOCH_BINARY]}
epoch=$(realpath "${2:-$here/../../../bin/epoch2d}")
top=$here/run_$label

# run DIR NPROCS: run EPOCH in DIR on input.deck
run() {
  ln -s "$(dirname "$epoch")/../src" "$1/src"
  (cd "$1" && echo . | mpirun -n "$2" "$epoch" > epoch.log 2>&1) \
    || { echo "   FAILED, see $1/epoch.log"; exit 1; }
  grep -h "Final runtime" "$1/epoch.log" || true
}

rm -rf "$top"
mkdir -p "$top/full"
cp "$here/input.deck" "$top/full/input.deck"
echo "== full [$label, 2 ranks]"
run "$top/full" 2
[ -f "$top/full/0005.sdf" ] || { echo "   no full/0005.sdf"; exit 1; }

for n in 2 1; do
  rdir=$top/restart$n
  mkdir -p "$rdir"
  sed 's/^begin:control$/begin:control\n  restart_snapshot = 5/' \
    "$here/input.deck" > "$rdir/input.deck"
  ln -s ../full/0005.sdf "$rdir/0005.sdf"
  echo "== restart from 0005 [$n ranks]"
  run "$rdir" "$n"
done
