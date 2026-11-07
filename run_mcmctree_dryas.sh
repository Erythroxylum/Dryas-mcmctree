#!/usr/bin/env bash
set -euo pipefail

# Run approximate likelihood preparation once and two MCMC chains for each of
# nine chromosome analyses created by setup_mcmctree_dryas.R.
# Usage: bash run_mcmctree_dryas.sh OUTPUT_DIR
# Optional: MCMCTREE=/path/to/mcmctree JOBS=4 bash run_mcmctree_dryas.sh OUTPUT_DIR

if [[ $# -ne 1 ]]; then
  echo "Usage: bash run_mcmctree_dryas.sh OUTPUT_DIR" >&2
  exit 2
fi

analysis_dir=$(cd "$1" && pwd)
mcmctree_bin=${MCMCTREE:-mcmctree}
jobs=${JOBS:-1}
force=${FORCE:-0}

command -v "$mcmctree_bin" >/dev/null 2>&1 || {
  echo "Cannot find mcmctree executable: $mcmctree_bin" >&2
  exit 1
}
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "JOBS must be a positive integer" >&2; exit 2; }
[[ "$force" =~ ^[01]$ ]] || { echo "FORCE must be 0 or 1" >&2; exit 2; }

run_prepare() {
  local chr_dir=$1
  if [[ -s "$chr_dir/prepare/out.BV" ]]; then
    echo "[$(basename "$chr_dir")] reusing existing out.BV"
    cp "$chr_dir/prepare/out.BV" "$chr_dir/chain1/in.BV"
    cp "$chr_dir/prepare/out.BV" "$chr_dir/chain2/in.BV"
    return 0
  fi
  echo "[$(basename "$chr_dir")] preparing approximate likelihood"
  (
    cd "$chr_dir/prepare"
    "$mcmctree_bin" mcmctree.ctl > mcmctree.console.log 2>&1 || {
      status=$?
      echo "MCMCTree preparation failed with status $status in $PWD" >&2
      tail -30 mcmctree.console.log >&2 || true
      exit "$status"
    }
    [[ -s out.BV ]] || { echo "Missing out.BV in $PWD" >&2; exit 1; }
    cp out.BV ../chain1/in.BV
    cp out.BV ../chain2/in.BV
    touch PREPARE_COMPLETE
  )
}

run_chain() {
  local chain_dir=$1
  if [[ "$force" != 1 && -s "$chain_dir/mcmc.txt" && -f "$chain_dir/CHAIN_COMPLETE" ]]; then
    echo "[$(basename "$(dirname "$chain_dir")")/$(basename "$chain_dir")] already complete"
    return 0
  fi
  echo "[$(basename "$(dirname "$chain_dir")")/$(basename "$chain_dir")] running"
  (
    cd "$chain_dir"
    rm -f mcmc.txt mcmctree.out FigTree.tre CHAIN_COMPLETE
    "$mcmctree_bin" mcmctree.ctl > mcmctree.console.log 2>&1 || {
      status=$?
      echo "MCMCTree chain failed with status $status in $PWD" >&2
      tail -30 mcmctree.console.log >&2 || true
      exit "$status"
    }
    [[ -s mcmc.txt ]] || { echo "Missing mcmc.txt in $PWD" >&2; exit 1; }
    touch CHAIN_COMPLETE
  )
}

export -f run_prepare run_chain
export mcmctree_bin force

# Bash arrays work on macOS without GNU sort -z or find depth options.
shopt -s nullglob
chromosomes=("$analysis_dir"/chr??)
[[ ${#chromosomes[@]} -eq 9 ]] || { echo "Expected nine chromosome folders" >&2; exit 1; }
for chr_dir in "${chromosomes[@]}"; do
  if ! grep -Eq '^usedata = 0$' "$chr_dir/chain1/mcmctree.ctl"; then
    run_prepare "$chr_dir"
  fi
done
printf '%s\0' "${chromosomes[@]/%//chain1}" "${chromosomes[@]/%//chain2}" |
  xargs -0 -n 1 -P "$jobs" bash -c 'run_chain "$1"' _
echo "All MCMCTree runs completed."
