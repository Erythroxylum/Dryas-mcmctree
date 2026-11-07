#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source config.sh
export REPRESENTATIVES CLADE_MEMBERSHIP P28_TREE
stage=${1:?Usage: bash workflow.sh STAGE DATASET SCENARIO}
export DATASET=${2:-both}
scenario=${3:-fossil}
if [[ "$DATASET" == both ]]; then
  bash workflow.sh "$stage" 15tip "$scenario"
  bash workflow.sh "$stage" p28 "$scenario"
  exit
fi
[[ "$DATASET" == 15tip || "$DATASET" == p28 ]] || exit 64
out="$RESULTS/$DATASET/$scenario"
[[ ${PRIOR_ONLY:-0} == 1 ]] && out="${out}_prior"
case "$stage" in
 setup)
  if [[ "$DATASET" == p28 ]]; then align=$ALIGNMENTS_P28; trees=.; else align=$ALIGNMENTS_15; trees=$TREES_15; fi
  Rscript setup_mcmctree_dryas.R "$align" "$trees" "$out" "$scenario" ;;
 run) bash run_mcmctree_dryas.sh "$out" ;;
 summarize) Rscript summarize_mcmctree_dryas.R "$out" ;;
 export) Rscript export_trees.R "$out" "$DATASET" "$REPRESENTATIVES" ;;
 mcc)
  for kind in full dryas; do
    "${TREEANNOTATOR:-treeannotator}" -burnin 0 -heights median "$out/exports/$kind.nex" "$out/exports/${kind}_MCC.tre"
  done ;;
 *) echo 'Stages: setup run summarize export mcc' >&2; exit 64 ;;
esac
