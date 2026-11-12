#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source config_majority.sh
export DATASET=p28 REPRESENTATIVES CLADE_MEMBERSHIP P28_TREE
stage=${1:?Usage: bash workflow_majority.sh prepare ROOTED_TREE_DIR | setup | run | summarize | export | mcc}
out="$RESULTS/p28/fossil"
[[ ${PRIOR_ONLY:-0} == 1 ]] && out="${out}_prior"
case "$stage" in
 prepare)
  [[ $# -eq 2 ]] || { echo "Supply rooted chromosome tree directory" >&2; exit 64; }
  Rscript prepare_majority_dating.R "$2" "$BASE_P28_TREE" "$REPRESENTATIVES" "$CLADE_MEMBERSHIP" "$MAJORITY_TOPOLOGY_DIR"
  ;;
 setup)
  [[ -d "$ALIGNMENTS_P28" ]] || { echo "Alignment directory not found: $ALIGNMENTS_P28. Override ALIGNMENTS_P28 with the existing filtered p28 alignment directory." >&2; exit 1; }
  [[ -s "$P28_TREE" ]] || { echo "Run prepare first; missing $P28_TREE" >&2; exit 1; }
  Rscript setup_mcmctree_dryas.R "$ALIGNMENTS_P28" "." "$out" fossil
  cp "$MAJORITY_TOPOLOGY_DIR/p28_majority_clade_audit.tsv" "$out/"
  cp "$MAJORITY_TOPOLOGY_DIR/p28_majority_resolution_choices.tsv" "$out/"
  ;;
 run) bash run_mcmctree_dryas.sh "$out" ;;
 summarize) Rscript summarize_mcmctree_dryas.R "$out" ;;
 export) Rscript export_trees.R "$out" p28 "$REPRESENTATIVES" ;;
 mcc)
  for kind in full dryas; do
   "${TREEANNOTATOR:-treeannotator}" -burnin 0 -heights median "$out/exports/$kind.nex" "$out/exports/${kind}_MCC.tre"
  done ;;
 *) echo "Stages: prepare setup run summarize export mcc" >&2; exit 64 ;;
esac
