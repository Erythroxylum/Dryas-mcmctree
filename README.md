# Unified PanDryas dating: 15-tip and p28

## Current status and next step

This repository runs independent chromosome dating analyses for the 15-tip and p28 datasets, exports 900 sampled trees per dataset, and documents summarization in the external TreeAnnotator app. For a new analysis follow Sections 1–3, then 5. For completed runs, start with diagnostics in Section 3; do not repeat setup. Matching prior-only analyses are described in Section 4.

The supplied local MCC files are `MCC_15tip_15.tre` (15 tips) and `MCC_28_32tip.tre` (32 tips). Both include four outgroups. They are single annotated summary trees, not the 900-tree posterior inputs; analysis-result trees are not included in this repository upload.

Commands below assume the `fossil` scenario and this project directory:

```bash
cd '/Users/dawsonwhite/Library/CloudStorage/OneDrive-TexasA&MAgriLife/Dryas/panDryas/pyrad/pandryas_mcmctree'
```

## 1. Software and inputs

Requires R packages `ape` and `coda`, Python 3 for optional input filtering, PAML `mcmctree`, and BEAST TreeAnnotator for MCC summaries.

```bash
Rscript -e 'install.packages(c("ape", "coda"), repos="https://cloud.r-project.org")'
```

The current workflow uses **separate alignment sources**:

- `chromosome_alignments_trimmed/`: nine original 15-tip chromosome PHYLIP alignments (11 Dryas and four outgroups).
- `chromosome_alignments_sitefiltered/`: nine already site-filtered p28 PHYLIP alignments (28 Dryas and four outgroups). These were copied into this project from the previous p28 analysis. The chromosome 1 header is `32 374743`.
- `iqtree_chromosomes_trimmed/`: nine original 15-tip IQ-TREE `.treefile` files. IQ-TREE result files are not sequence alignments.
- `pop28_fixed_samples.tre`: rooted, bifurcating 32-tip p28 topology.
- `pop28_representatives.tsv`: 28 population-to-sample mappings.
- `clade_membership.tsv`: population metadata and explicit calibration membership.

Use these settings in `config.sh`:

```bash
ALIGNMENTS_15="chromosome_alignments_trimmed"
ALIGNMENTS_P28="chromosome_alignments_sitefiltered"
TREES_15="iqtree_chromosomes_trimmed"
P28_TREE="pop28_fixed_samples.tre"
REPRESENTATIVES="pop28_representatives.tsv"
CLADE_MEMBERSHIP="clade_membership.tsv"
RESULTS="results_unified"
```

The current clade table has five tab-separated columns:

```text
sampleID    population    tip.label    clade    calibration_clade
```

`clade` retains detailed lineage names. **Setup must read `calibration_clade`** for the three calibration groups: `drummondii`, `AIH`, and `other_core`. Merely adding that column does not update an old setup script that still reads `clade`.

Current p28 membership is nine AIH samples (two ajanensis, two alaskensis, two hookeriana, three integrifolia), one drummondii sample, and 18 other-core samples. Setup checks monophyly. Membership is explicit and must not be inferred from sample prefixes. `tip.label` stores geographic display labels; its presence alone does not make the exporter or plotting script use it.

**Do not run `prepare_inputs.R` or the p28 filter again for the current inputs.** The original preparation route assumed alignments containing all required representatives; the 15-tip files cannot supply the additional p28 samples. Renaming samples cannot reconstruct missing sequences. Both datasets retain their own representatives and sites; the 15-tip analysis is not obtained by pruning p28.

## 2. Setup and run — for a new analysis only

For fresh output directories:

```bash
bash workflow.sh setup both fossil
JOBS=4 bash workflow.sh run both fossil
```

`both` processes 15tip and p28 in sequence. Each dataset has nine chromosomes and two chains per chromosome. `JOBS` controls runner concurrency. Sampling settings are 100000 burn-in iterations, sampling every 20 iterations, and 100000 retained samples per chain.

Outputs are under `results_unified/15tip/fossil/` and `results_unified/p28/fossil/`, with chain files at `chr01/chain1/mcmc.txt`, etc. The message `Crown-Dryas calibration: none` is expected for the fossil scenario: root, Cercocarpus, core-Dryas and AIH constraints still apply.

Setup requires an empty output directory. The runner skips chains marked complete; this is not continuation of a chain from its last sample. Do not change inputs or settings and reuse old completion flags. For changed analyses, use a new `RESULTS` directory in `config.sh` and keep previous results intact.

## 3. Check convergence now

### Generate summaries

```bash
bash workflow.sh summarize both fossil
```

Each dataset receives:

| File | Contents |
|---|---|
| `mcmctree_convergence_summary.tsv` | Per-chromosome minimum ESS and maximum classical R-hat |
| `mcmctree_parameter_convergence.tsv` | ESS for each chain and R-hat for every parameter |
| `mcmctree_chain_posterior_summary.tsv` | Separate chain estimates and equal-tailed 95% intervals |
| `mcmctree_pooled_posterior_summary.tsv` | Estimates pooling the two chains within each chromosome |
| `mcmctree_posterior_summary.tsv` | Compatibility copy of the pooled summary |

Compare the two chains **within each chromosome and dataset**. Do not calculate convergence across different chromosomes or combine their chains to improve ESS.

Use finite classical R-hat <= 1.05 and ESS >= 200 in **each** chain as initial screening criteria. These are screening thresholds, not proof of convergence. Inspect all sampled parameters, including node ages, likelihood where present, and clock/rate parameters. NA or infinite diagnostics require investigation; they are not passes. Constant parameters can produce undefined diagnostics and should be checked explicitly.

Print all parameter failures with this command:

```bash
Rscript - <<'RS'
for (dataset in c("15tip", "p28")) {
  root <- file.path("results_unified", dataset, "fossil")
  cat("\nDATASET:", dataset, "\n")
  print(read.delim(file.path(root, "mcmctree_convergence_summary.tsv")), row.names=FALSE)
  d <- read.delim(file.path(root, "mcmctree_parameter_convergence.tsv"))
  pass <- with(d, is.finite(Rhat) & Rhat <= 1.05 &
    is.finite(chain1_ESS) & chain1_ESS >= 200 &
    is.finite(chain2_ESS) & chain2_ESS >= 200)
  cat("Parameters failing the initial screen:", sum(!pass), "of", nrow(d), "\n")
  if (any(!pass)) print(d[!pass, ], row.names=FALSE)
}
RS
```

### Check sample counts and inspect traces

A count of 100001 rows was observed in a completed run. The code reports unexpected counts as warnings rather than rejecting the plot. Inspect the iteration column and control settings to explain the difference; do not silently delete a row.

The following command reports retained sample counts, checks matching parameter columns and finite numeric values, and creates an overlay-trace PDF for each dataset. It plots every numeric parameter except the iteration counter. Only the displayed traces are thinned to at most 5000 points per chain; summary diagnostics use all retained draws.

```bash
Rscript - <<'RS'
read_chain <- function(path) {
  x <- read.table(path, header=TRUE, check.names=FALSE, comment.char="")
  if (!nrow(x)) stop("Empty chain: ", path)
  if (nrow(x) != 100000L) warning("Check retained sample count: ", path, " (", nrow(x), ")")
  x <- x[, vapply(x, is.numeric, logical(1)), drop=FALSE]
  x <- x[, !tolower(names(x)) %in% c("gen", "iter", "iteration", "generation", "sample", "index"), drop=FALSE]
  if (!ncol(x) || any(!is.finite(as.matrix(x)))) stop("Missing/nonfinite numeric parameters: ", path)
  x
}
for (dataset in c("15tip", "p28")) {
  root <- file.path("results_unified", dataset, "fossil")
  pdf(file.path(root, "mcmctree_trace_overlays.pdf"), width=11, height=8.5)
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  for (chr in 1:9) {
    paths <- file.path(root, sprintf("chr%02d", chr), c("chain1", "chain2"), "mcmc.txt")
    x <- lapply(paths, read_chain)
    if (!identical(names(x[[1]]), names(x[[2]]))) stop("Parameter mismatch: ", dataset, " chr", chr)
    cat(dataset, sprintf("chr%02d", chr), "retained rows:", nrow(x[[1]]), nrow(x[[2]]), "\n")
    ii <- unique(round(seq(1, nrow(x[[1]]), length.out=min(5000, nrow(x[[1]])))))
    for (p in names(x[[1]])) {
      yr <- range(x[[1]][[p]], x[[2]][[p]])
      plot(ii, x[[1]][[p]][ii], type="l", col=adjustcolor("#0072B2", .55),
           ylim=yr, xlab="Retained draw", ylab=p,
           main=paste(dataset, sprintf("chr%02d", chr), p))
      lines(ii, x[[2]][[p]][ii], col=adjustcolor("#D55E00", .55))
      legend("topright", c("chain1", "chain2"), col=c("#0072B2", "#D55E00"), lty=1, bty="n", cex=.7)
    }
  }
  dev.off()
}
RS

open results_unified/15tip/fossil/mcmctree_trace_overlays.pdf
open results_unified/p28/fossil/mcmctree_trace_overlays.pdf
```

Look for stable, overlapping chains with no sustained drift, persistent separation, or long stuck periods. Compare chain-specific medians and intervals where traces or diagnostics are concerning. MCMCTree's configured burn-in precedes the retained `mcmc.txt` samples; these commands do not discard another fraction automatically.

If checks fail, diagnose the affected chromosome/parameter before export. Longer independent runs, improved tuning, or a revised burn-in may be needed. Increasing thinning or simply discarding samples is not a general convergence fix. If any additional burn-in is justified, apply it consistently to diagnostics and export; the current commands do not implement that change.

## 4. Prior-only runs

Run the matching prior scenario to inspect the effective age prior induced jointly by calibrations and the tree prior. Prior-only output also needs convergence assessment. These are additional analyses, not a substitute for posterior diagnostics.

```bash
PRIOR_ONLY=1 bash workflow.sh setup both fossil
PRIOR_ONLY=1 JOBS=4 bash workflow.sh run both fossil
PRIOR_ONLY=1 bash workflow.sh summarize both fossil
```

Outputs use `fossil_prior` instead of `fossil`. To reuse the diagnostic and trace commands above, change the final path component from `"fossil"` to `"fossil_prior"`. Do not repeat setup if these directories already contain an analysis. Compare corresponding node-age prior and posterior distributions within chromosomes to assess how much the data inform ages.

## 5. Export and summarize accepted posterior runs

After convergence review, generate the posterior tree files. This command runs only R; it does not invoke TreeAnnotator:

```bash
bash workflow.sh export both fossil
```

### TreeAnnotator app (primary route)

Open TreeAnnotator separately. Run it once for each input you need:

| Input Tree File | Suggested Output File |
|---|---|
| `results_unified/15tip/fossil/exports/full.nex` | `results_unified/15tip/fossil/exports/full_MCC.tre` |
| `results_unified/15tip/fossil/exports/dryas.nex` | `results_unified/15tip/fossil/exports/dryas_MCC.tre` |
| `results_unified/p28/fossil/exports/full.nex` | `results_unified/p28/fossil/exports/full_MCC.tre` |
| `results_unified/p28/fossil/exports/dryas.nex` | `results_unified/p28/fossil/exports/dryas_MCC.tre` |

Select **burn-in 0**, **Maximum clade credibility tree**, and **Median node heights**. Use posterior probability limit 0 to retain annotations for all nodes. Each input contains 900 trees (100 per chromosome, 50 per chain). Select the `.nex` file as input, not an existing MCC file. An app bundle need not be a command-line executable.

On macOS, open the input folders in Finder:

```bash
open results_unified/15tip/fossil/exports
open results_unified/p28/fossil/exports
```

### Optional command-line TreeAnnotator

If a working command-line executable is installed, the existing workflow stage performs the same summaries:

```bash
command -v treeannotator
bash workflow.sh mcc both fossil
```

Or supply the actual executable path (not merely the `.app` directory):

```bash
TREEANNOTATOR='/absolute/path/to/treeannotator' bash workflow.sh mcc both fossil
```

Skip this stage when using the app manually. Rerunning MCMCTree is unnecessary for either route.

After export and TreeAnnotator summarization, each dataset's `exports/` directory contains `full.nex`, `dryas.nex`, `full_MCC.tre`, and `dryas_MCC.tre`. The full trees have 15 or 32 tips; Dryas-only trees have 11 or 28 tips. MCMCTree estimates ages on the supplied topology; it does not infer a posterior distribution of topologies.

The exporter selects 50 retained draws per chain per chromosome: 900 total, equally weighted across nine chromosomes. This is an **equal-weight chromosome mixture**, not a joint multilocus posterior. Convergence must be assessed before this small visualization sample is exported.

The exporter maps sampled ages using PAML's numbered topology in each `mcmctree.out`, rather than assuming PAML and ape node numbers agree. It stops if the required topology format cannot be recovered. If that occurs, inspect the output format and update the parser; do not bypass the check.

Node-age values in raw samples and summary TSVs are in units of 100 Ma: multiply age values by 100 for Ma. Do not multiply likelihoods or rate parameters by 100. Exported tree branch lengths are already in Ma.

TreeAnnotator uses `-burnin 0 -heights median` because the exported samples are already post-burn-in. Preserve its annotations, including `height_95%_HPD`. These HPD intervals differ from the equal-tailed intervals in the per-chromosome summary tables. The mixture's intervals are also distinct from individual-chromosome intervals.

## 6. Figure preparation

Run the population tree/map plot from the project root. Supply either the full 32-tip annotated MCC or the pruned 28-tip Dryas MCC:

```bash
Rscript -e 'install.packages(c("ape", "sf", "maps", "BiocManager"), repos="https://cloud.r-project.org"); BiocManager::install("treeio", ask=FALSE)'
Rscript plot_p28_equalearth_hpd.R MCC_28_32tip.tre
```

The second command accepts an explicit tree path; for a newly summarized fossil analysis use `results_unified/p28/fossil/exports/full_MCC.tre`.

Place `Dryas_sampledata_popgen.csv` in the project root, or supply its path as the second argument:

```bash
Rscript plot_p28_equalearth_hpd.R /path/to/MCC_28_28tip.tre /path/to/Dryas_sampledata_popgen.csv
```
 Required columns are sampleID, s170-p28, s170_BPP, s380, Longitude, Latitude, and Epithet. If labels or species assignments differ, provide `population_species_key.csv` with tree_tip,population,species, one row per population.

The script reads node HPDs directly from the selected MCC, prunes four named outgroups, matches intervals by descendant populations, and scales the age axis to include every HPD. Plotted ages come from the supplied branch lengths, supporting mean or median summaries. Branch lengths are already in Ma. No separate old HPD CSV is required.

Outputs are `figures/figure3_p28_species_map.pdf`, `.png`, and node-age, color, and sample-audit CSVs. This produces the dated population tree/map component, not the BPP model-diagram panel.

This revision has not been executed in the editing environment because R is unavailable. Run the command locally before treating the outputs as validated.

## 7. Scenarios and calibrations

Use the scenario name as the last workflow argument:

| Scenario | Additional Dryas constraints |
|---|---|
| `fossil` | Core-Dryas minimum 6.0 Ma; crown-AIH minimum 3.4 Ma; no BPP crown constraint |
| `combined` | Fossil constraints plus chromosome-specific BPP crown-Dryas interval |
| `baseline` | BPP crown-Dryas interval; neither new fossil constraint |
| `japan_stem` | Japanese minimum moved to crown Dryas (the core stem divergence); AIH minimum retained; no BPP crown constraint |
| `outgroup_only` | Neither new fossil constraint nor BPP crown constraint |

All retain root B(34.07,48.70 Ma), stem-Cercocarpus L(34.07 Ma), clock=2, and the setup script's model/rate priors. Calibration parameters use 100 Ma units: 6.0 Ma is 0.060, 3.4 Ma is 0.034. PAML L() and B() use default soft-bound parameters; they are not hard truncations.

Japanese fossil placement at core crown and Beaver Pond placement at crown AIH are explicit morphological/phylogenetic assumptions. These alternative scenarios evaluate some calibration sensitivity; a minimum fossil age alone does not establish crown membership.

Sources:

- Japanese fossil: https://doi.org/10.2517/prpsj.250006
- Arctic macroflora: https://doi.org/10.26879/1121
- Beaver Pond age context: https://doi.org/10.5194/cp-15-1063-2019

Use separate output directories for each scenario, and matching prior-only runs when comparing effective priors.

## 8. Repository contents and large files

The repository contains the setup, runner, summarizer, exporter, workflow wrapper, configuration, optional preparation/filtering scripts, but no analysis-result trees. Alignments, input topologies, and local population/membership tables must be supplied separately. The legacy preparation and filtering scripts are retained for alternative input routes; they are not required with the current separate alignment directories. Final geographic plotting remains a separate workflow.

`.gitignore` excludes `mcmc.txt` and `mcmctree.out` at every directory depth, plus alignment directories, large approximate-likelihood files, logs, and completion flags. Small summary tables and trees remain eligible for version control. Review `git status` before committing generated results.

**Keep the ignored chain files locally or in a separate archive.** The exporter needs both `mcmc.txt` and `mcmctree.out`; ignoring them is not an instruction to delete them. A fresh clone cannot recreate posterior exports without archived chain outputs or rerunning the analyses. Ignore rules do not remove files already tracked by Git.

## 9. Validation status

The runs were reported complete and parameters well mixed on the user's Mac; two TreeAnnotator MCC files were supplied afterward. The repository version incorporates the documented alignment-path and calibration-column fixes. Shell syntax and Python syntax are checked before upload. R/PAML execution of this packaged revision has not been tested here; retain the convergence summaries and inspect new results after any changed settings.


## PhyParts chromosome concordance pies

Generate the pie table locally from the PhyParts output prefix (requires .hist, .concon.tre, and .node.key):

```bash
Rscript prepare_phyparts_pies.R "$BASE/phyparts_p28"
Rscript plot_p28_equalearth_hpd.R results_unified/p28/fossil/exports/MCC_28_28tip.tre Dryas_sampledata_popgen.csv
```

The converter creates data/phyparts_p28_pies.csv. Derived counts are supplied locally, not included in the repository. A custom CSV path can be supplied as the converter's second argument and the plotter's third argument.

Black: concordant/9. White: largest alternative count/9. Gray: (total conflicting trees - largest alternative count)/9. Alternative counts are not summed. Counts come from PhyParts outputs and must be nonnegative and account for all nine chromosomes. Uninformative signal causes a stop, because it would require a fourth category.

Clades are joined by descendant-tip sets and the population rename key. The root and the 27-tip non-drum clade forced by drum rooting are not given pies. Ages and HPDs are read from the MCC as before. Plot outputs include a _concordance.csv audit with node mappings. R plotting execution is not verified in the editing environment, where R is unavailable.


## Fossil leaf markers in Figure 3

The plot uses the supplied assets/Fig1-leaf.svg at the crown of non-drummondii Dryas (6.0 Ma minimum) and crown AIH (3.4 Ma minimum). These markers assume the fossil calibration scenario. The script checks explicit clade membership, writes a _fossil_calibrations.csv audit, and offsets each leaf from the node and concordance pie. Adjust LEAF_HEIGHT, LEAF_DX, and LEAF_DY near the top of the plot script for placement.

Install the additional SVG renderer in the active Conda environment:

```bash
conda install -c conda-forge r-rsvg
```

The original SVG is retained; rsvg renders a high-resolution transparent bitmap for embedding in base-R PDF/PNG output. The concordance annotation is omitted. Placement has not been rendered in R in the editing environment.


## All chromosome and majority-rule topology figures

From the repository root, supply the directory containing the nine individual CASTER trees rooted on drum:

```bash
Rscript plot_p28_topology_summaries.R "$BASE/phyparts_input_rooted"
```

Requires only ape. The directory must contain one tree per chromosome (chr1.tre or s170-p28.ch1.drum_rooted.tre naming is accepted). All trees must have the same 28 tips and drum as one of two root daughters.

Outputs in figures/topology_summaries/:
- p28_chromosome1 through p28_chromosome9: each has .pdf and .png figures, a .nwk topology, and a _pies.csv table.
- p28_all_chromosomes.pdf: nine pages, one chromosome per page, in chromosome order.
- p28_majority_rule.pdf and .png, .nwk topology, and _pies.csv.

The summary is a rooted majority-rule consensus: clades must occur in at least five of the nine trees; unresolved relationships remain polytomies. All chromosomes receive equal weight. These are cladograms with arbitrary display lengths, no date axis, and no fossil markers.

Pie counts are recalculated for each reference. Black counts exact matching rooted clades; white counts the most frequent incompatible clade; gray is the number of distinct conflicting chromosome trees minus that maximum. Rooted clades conflict if they overlap and neither contains the other. This direct method is used only for the complete, identical taxon sets supplied here. Its complete alternative histograms reproduced all 26 clades of the supplied PhyParts primary-reference results. Uninformative signal triggers an error. Root and forced non-drum sister-clade pies are omitted. Pie colors and omission of the concordance legend match Figure 3.

To generate the existing dated figure and both companions together, provide the chromosome directory as the fourth argument:

```bash
Rscript plot_p28_equalearth_hpd.R results_unified/p28/fossil/exports/MCC_28_28tip.tre Dryas_sampledata_popgen.csv data/phyparts_p28_pies.csv "$BASE/phyparts_input_rooted"
```

The R plotting code has not been executed in the editing environment because R is unavailable.


### Fully bifurcating MCC summary

The same topology-summary command also exports p28_mcc.pdf, .png, .nwk, and _pies.csv, plus p28_mcc_scores.csv for all nine candidates. p28_summaries.pdf contains the majority-rule consensus and MCC on two separate pages.

MCC selection maximizes the sum of log clade frequencies across the candidate's internal clades, using frequencies among the nine chromosomes and excluding the whole-tree root. The forced non-drum clade has frequency 1 and does not affect the score. All candidates must be fully bifurcating after redundant single-child nodes are collapsed. Exact score ties within 1e-10 are reported; the lowest chromosome number is selected deterministically.

This is a chromosome-frequency MCC, not a posterior-probability tree: the inputs are nine chromosome point estimates. It selects one observed topology and does not arbitrarily resolve majority-rule polytomies. The majority-rule output is retained. Neither summary receives calendar dates or fossil markers. Independently calculated scores for the supplied nine trees select chromosome 1; R rendering remains unverified in the editing environment.


## Separate dating on a resolved majority-rule topology

This p28-only comparison has its own config_majority.sh and workflow_majority.sh. It writes to results_majority_resolved/p28/fossil/, preserving the original results_unified runs. It uses fossil calibrations only, with no BPP crown-Dryas constraint. Two chains and all nine chromosome alignments retain the original clock, rate priors, and sampling settings.

MCMCTree requires a binary reference. prepare_majority_dating.R retains every clade occurring in at least five chromosomes, then adds observed compatible clades in descending frequency order. Equal-frequency ties use lexical clade keys (deterministic, not evidence for preferring the winner). The exact added clades and tied candidates are recorded in topologies/majority_resolved/p28_majority_resolution_choices.tsv. For the supplied nine trees, the three resolutions have 4/9, 3/9, and 1/9 support. This is a resolved majority-rule reference, not the raw consensus or an MCC tree.

The 28 population tips are mapped back to sample IDs with pop28_representatives.tsv. The existing pop28_fixed_samples.tre supplies the four-outgroup scaffold only; its Dryas subtree is replaced. Exact core-Dryas and AIH calibration membership is verified against clade_membership.tsv. The 6.0 Ma core-Dryas and 3.4 Ma AIH minima and original outgroup/root constraints are retained.

From the repository root:

```bash
# BASE is the CASTER chromosome-results directory used for PhyParts.
bash workflow_majority.sh prepare "$BASE/phyparts_input_rooted"

# Review the three resolutions before launching the dating analysis.
cat topologies/majority_resolved/p28_majority_resolution_choices.tsv

bash workflow_majority.sh setup
JOBS=4 bash workflow_majority.sh run
bash workflow_majority.sh summarize
```

The default filtered alignment folder is p28_chromosome_alignments_sitefiltered, matching the supplied repository inventory. If yours is named chromosome_alignments_sitefiltered instead, use:

```bash
ALIGNMENTS_P28=chromosome_alignments_sitefiltered bash workflow_majority.sh setup
```

Setup requires a new empty results directory. Prepare also stops if its generated topology files already exist. Do not delete completed runs or repeat setup in them. Run can be invoked again to skip completed chains. The new topology needs fresh likelihood derivatives and branch-length optimization; no old chromosome-4 in.BV/out.BV files are copied.

Inspect convergence and trace plots using the original README guidance with results_majority_resolved replacing results_unified. After accepting diagnostics:

```bash
bash workflow_majority.sh export
open results_majority_resolved/p28/fossil/exports
```

Run TreeAnnotator separately on full.nex and dryas.nex (burn-in 0, MCC topology, median heights) or use the optional CLI stage:

```bash
bash workflow_majority.sh mcc
```

Optional matching prior-only runs use PRIOR_ONLY=1 for setup, run, summarize, and export, writing fossil_prior. Prepare is shared and runs only once. Change MAJORITY_TOPOLOGY_DIR and RESULTS_MAJORITY to new directories to explore different resolution policies; never overwrite the previous analysis.

Validation in the editing environment: shell syntax checked; the equivalent clade-selection calculation on all nine supplied CASTER trees preserved all majority clades, resolved the three remaining branches using observed compatible clades, and retained AIH. R/PAML execution, local sample mappings, and outgroup grafting must be checked by the workflow on the user's files; no new dating run has been executed here.
