# Separate p28 fossil-only dating on a resolved majority-rule reference.
# Override paths through environment variables when needed.
ALIGNMENTS_P28="${ALIGNMENTS_P28:-p28_chromosome_alignments_sitefiltered}"
REPRESENTATIVES="${REPRESENTATIVES:-pop28_representatives.tsv}"
CLADE_MEMBERSHIP="${CLADE_MEMBERSHIP:-clade_membership.tsv}"
BASE_P28_TREE="${BASE_P28_TREE:-pop28_fixed_samples.tre}"
MAJORITY_TOPOLOGY_DIR="${MAJORITY_TOPOLOGY_DIR:-topologies/majority_resolved}"
P28_TREE="$MAJORITY_TOPOLOGY_DIR/p28_majority_resolved.samples32.tre"
RESULTS="${RESULTS_MAJORITY:-results_majority_resolved}"
