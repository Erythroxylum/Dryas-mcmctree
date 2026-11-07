# Run from this folder in RStudio: source("plot_p28_equalearth_hpd.R")
# Uses your pruned all-chromosome MCC chronogram without modifying topology or dates.
packages <- c("ape", "sf", "maps", "treeio")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) stop("Install required packages first: ", paste(missing, collapse=", "))
library(ape)
library(sf)

args <- commandArgs(trailingOnly=TRUE)
TREE_FILE <- if (length(args)) args[1] else "MCC_28_32tip.tre"
METADATA_FILE <- if (length(args) >= 2L) args[2] else "Dryas_sampledata_popgen.csv"
OUTPUT_PREFIX <- "figures/figure3_p28_species_map"
dir.create("figures", showWarnings=FALSE)
AGE_MULTIPLIER <- 1 # Supplied exported MCC branch lengths are already in Ma.
CENTRAL_LONGITUDE <- -20 # Standard seam; no coastline shift.
MIN_LATITUDE <- 27
MAX_AGE <- 26
JITTER_KM <- 35 # Maximum projected displacement, for display only.
JITTER_SEED <- 42
MAP_SCALE <- 1.12 # Enlarge the map by 12 percent.

POP_COLUMN <- "s170-p28"
SPECIES_COLUMN <- "Epithet" # Case-insensitive column lookup.
# Optional CSV: tree_tip,population,species. Useful if tree names are older.
KEY_FILE <- "population_species_key.csv"
CONNECTION_ALPHA <- 0.17
TREE_HEIGHT <- 0.65 # In map-radius units.
TIP_BASELINE <- 0.65 # Move down to reduce the gap; allow room for tip labels.
REORDER_FOR_GEOGRAPHY <- TRUE # Rotates sister branches only; preserves topology.

if (!file.exists(TREE_FILE)) stop("Tree file not found: ", TREE_FILE)
if (!file.exists(METADATA_FILE)) stop("Metadata file not found: ", METADATA_FILE, "\nSupply its path as the second argument: Rscript plot_p28_equalearth_hpd.R TREE METADATA")
annotated <- treeio::read.beast(TREE_FILE)
full_tree <- ape::as.phylo(annotated)
annotations <- as.data.frame(annotated@data)
# treeio may store node IDs as character strings. Numeric IDs are required
# for internal-node filtering and positional indexing into the phylogeny.
annotations$node <- suppressWarnings(as.integer(as.character(annotations$node)))
if (anyNA(annotations$node) || anyDuplicated(annotations$node) ||
    any(!annotations$node %in% seq_len(Ntip(full_tree) + full_tree$Nnode)))
  stop("MCC annotations contain missing, duplicate, or invalid node IDs.")
full_n <- Ntip(full_tree)
full_children <- split(full_tree$edge[,2], full_tree$edge[,1])
full_desc <- function(node) {
  if (node <= full_n) return(full_tree$tip.label[node])
  unlist(lapply(full_children[[as.character(node)]], full_desc), use.names=FALSE)
}
outgroups <- c("Cermon_DW847", "Chafol_1611", "Pursta_DW848", "Purtri_DW846b")
present_outgroups <- intersect(outgroups, full_tree$tip.label)
if (!((full_n == 32L && length(present_outgroups) == 4L) ||
      (full_n == 28L && length(present_outgroups) == 0L)))
  stop("Supply either the 32-tip full MCC or the 28-tip Dryas-only MCC.")
keep <- vapply(annotations$node, function(z) z > full_n &&
  !any(full_desc(z) %in% outgroups), logical(1))
ann <- annotations[keep, , drop=FALSE]
intervals <- ann[["height_0.95_HPD"]]
if (is.null(intervals)) intervals <- ann[["height_95%_HPD"]]
if (is.null(intervals)) stop("Missing node height HPD annotations in MCC.")
hpd_source <- data.frame(
  descendants=vapply(ann$node, function(z) paste(sort(full_desc(z)),collapse=";"), ""),
  mean_Ma=max(ape::node.depth.edgelength(full_tree)[seq_len(full_n)]) - ape::node.depth.edgelength(full_tree)[ann$node],
  lower_Ma=vapply(intervals, function(z) as.numeric(z)[1], 0.0),
  upper_Ma=vapply(intervals, function(z) as.numeric(z)[2], 0.0))
MAX_AGE <- ceiling(max(hpd_source$upper_Ma)/2)*2
tr <- if (length(present_outgroups)) drop.tip(full_tree, present_outgroups) else full_tree
if (Ntip(tr) != 28L) stop("Expected the 28-tip Dryas-only tree. Prune outgroups explicitly first.")
if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0))
  stop("The input must contain finite, nonnegative dated branch lengths.")
tr$edge.length <- tr$edge.length * AGE_MULTIPLIER
if (!is.rooted(tr) || !is.ultrametric(tr, tol = 1e-6)) stop("Input tree must be rooted and ultrametric.")

d <- read.csv(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
required <- c("sampleID", POP_COLUMN, "s170_BPP", "s380", "Longitude", "Latitude")
if (!all(required %in% names(d))) stop("Missing metadata columns: ", paste(setdiff(required, names(d)), collapse=", "))
d$in_s170 <- toupper(trimws(as.character(d$s170_BPP))) %in% c("TRUE", "T", "1", "YES")
d$in_s380 <- !is.na(d$s380) & nzchar(trimws(d$s380)) &
  !(toupper(trimws(d$s380)) %in% c("FALSE", "0", "NO"))
if (any(d$in_s170 & !d$in_s380)) stop("S170 includes samples not flagged in S380; inspect metadata.")
d <- d[d$in_s380, ]
if (anyDuplicated(d$sampleID)) stop("Duplicate sample IDs in S380.")
d$population <- trimws(d[[POP_COLUMN]])
d$Longitude <- as.numeric(d$Longitude); d$Latitude <- as.numeric(d$Latitude)
if (any(!is.finite(d$Longitude) | !is.finite(d$Latitude))) stop("Missing/invalid coordinates.")
if (any(abs(d$Longitude) > 180 | abs(d$Latitude) > 90)) stop("Coordinates outside longitude/latitude bounds.")
if (any(d$Latitude < MIN_LATITUDE)) stop("Some samples are south of the map boundary; inspect before excluding.")
if (any(is.na(d$population[d$in_s170]))) stop("Missing S170 population assignments.")
species_col <- names(d)[tolower(names(d)) == tolower(SPECIES_COLUMN)]
if (length(species_col) != 1) stop("Cannot uniquely find the epithet column.")
d$species <- trimws(as.character(d[[species_col]]))
# No automatic stripping of cf., or inference from old population prefixes.
key <- if (file.exists(KEY_FILE)) read.csv(KEY_FILE, stringsAsFactors=FALSE) else NULL
if (!is.null(key)) {
  stopifnot(all(c("tree_tip", "population", "species") %in% names(key)))
  if (anyDuplicated(key$tree_tip) || anyDuplicated(key$population) ||
      !setequal(key$tree_tip,tr$tip.label)) stop("Key must map all 28 tree tips one-to-one to populations.")
  tr$tip.label <- key$population[match(tr$tip.label,key$tree_tip)]
  # Map metadata labels using the same explicit key as the tree.
  ix <- match(d$population, key$tree_tip)
  mapped <- !is.na(ix)
  ambiguous <- mapped & d$population %in% key$population &
    d$population != key$population[ix]
  if (any(ambiguous, na.rm=TRUE))
    stop("Ambiguous key: a label is both an original and a different revised population.")
  d$population[mapped] <- key$population[ix[mapped]]
}
if (!setequal(unique(d$population[d$in_s170]), tr$tip.label)) {
  metadata_pops <- unique(d$population[d$in_s170])
  stop("Population mismatch after applying key.\nOnly in S170 metadata: ",
       paste(setdiff(metadata_pops, tr$tip.label), collapse=", "),
       "\nOnly in tree: ", paste(setdiff(tr$tip.label, metadata_pops), collapse=", "))
}
pop_species <- setNames(rep(NA_character_,length(tr$tip.label)),tr$tip.label)
for (pop in tr$tip.label) {
  values <- unique(na.omit(d$species[d$in_s170 & d$population == pop]))
  values <- values[nzchar(values)]
  if (!is.null(key)) {
    override <- key$species[match(pop,key$population)]
    if (!is.na(override) && nzchar(override)) values <- override
  }
  if (length(values) != 1) stop("Population ",pop," has ambiguous epithets: ",paste(values,collapse=", "),
    ". Resolve in your updated datasheet or the explicit population key.")
  pop_species[pop] <- values
}

# Shift longitudes so the 160 E seam becomes the usual +/-180 seam.
# Project the shifted coordinates at lon_0=0: equivalent to lon_0=-20.
crs_map <- "+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs"
map_points <- d
map_points$Longitude <- ((d$Longitude - CENTRAL_LONGITUDE + 180) %% 360) - 180
xy <- st_coordinates(st_transform(st_as_sf(map_points, coords=c("Longitude", "Latitude"), crs=4326), crs_map))
# Densified boundary of the complete northern cap; its sides narrow northward.
eps <- 1e-7
boundary_ll <- rbind(
  cbind(seq(-180+eps,180-eps,length.out=721),MIN_LATITUDE),
  cbind(180-eps,seq(MIN_LATITUDE,90,length.out=361)),
  cbind(seq(180-eps,-180+eps,length.out=721),90),
  cbind(-180+eps,seq(90,MIN_LATITUDE,length.out=361)))
boundary_ll <- rbind(boundary_ll,boundary_ll[1,])
boundary <- st_transform(st_sfc(st_polygon(list(boundary_ll)),crs=4326),crs_map)
box <- st_bbox(boundary)
radius <- max(abs(box[c("xmin","xmax")])) / MAP_SCALE
ybottom <- unname(box["ymin"])
map_top <- unname((box["ymax"]-ybottom)/radius)
#TIP_BASELINE <- map_top + 0.18
outline <- (st_geometry(boundary)-c(0,ybottom))/radius
set.seed(JITTER_SEED)
# Stable per-sample offsets independent of metadata row order.
o <- order(d$sampleID); theta <- runif(nrow(d),0,2*pi)
rj <- sqrt(runif(nrow(d))) * JITTER_KM * 1000
jx <- jy <- numeric(nrow(d)); jx[o] <- rj*cos(theta); jy[o] <- rj*sin(theta)
d$x_true <- xy[,1]/radius; d$y_true <- (xy[,2]-ybottom)/radius
d$x <- pmax(-MAP_SCALE,pmin(MAP_SCALE,(xy[,1]+jx)/radius))
d$y <- pmax(0,pmin(map_top,(xy[,2]+jy-ybottom)/radius))
# Keep jittered symbols inside the curved map boundary.
candidate <- st_as_sf(data.frame(x=d$x*radius,y=d$y*radius+ybottom),coords=c("x","y"),crs=crs_map)
outside <- lengths(st_intersects(candidate,boundary)) == 0
# Shrink offsets near the seam/boundary until they fit.
for (k in seq_len(12)) {
  if (!any(outside)) break
  d$x[outside] <- (d$x[outside]+d$x_true[outside])/2
  d$y[outside] <- (d$y[outside]+d$y_true[outside])/2
  candidate <- st_as_sf(data.frame(x=d$x*radius,y=d$y*radius+ybottom),coords=c("x","y"),crs=crs_map)
  outside <- lengths(st_intersects(candidate,boundary)) == 0
}
d$x[outside] <- d$x_true[outside]; d$y[outside] <- d$y_true[outside]

# Validate and record dates without moving nodes or refitting branch lengths.
n <- Ntip(tr); root <- setdiff(tr$edge[,1], tr$edge[,2])
children <- split(tr$edge[,2], tr$edge[,1])
descendants <- function(node) {
  if (node <= n) return(tr$tip.label[node])
  unlist(lapply(children[[as.character(node)]], descendants), use.names=FALSE)
}
depth <- node.depth.edgelength(tr)
root_age <- max(depth[seq_len(n)])
if (root_age > MAX_AGE) stop("Tree root exceeds plotting limits.")
ages <- root_age - depth
ages[seq_len(n)] <- 0
node_table <- data.frame(node=(n+1):(n+tr$Nnode), age_Ma=ages[(n+1):(n+tr$Nnode)],
  descendants=vapply((n+1):(n+tr$Nnode), function(z) paste(sort(descendants(z)), collapse=";"), ""))
# Intervals extracted from height_95%_HPD in the supplied all-chromosome MCC.
# Match by descendant populations, never by unstable node numbers.
hpd <- hpd_source
hpd$descendants <- vapply(strsplit(hpd$descendants,";",fixed=TRUE),function(labels) {
  if (!is.null(key)) {
    ix <- match(labels,key$tree_tip)
    if (anyNA(ix)) stop("HPD labels do not match the population key.")
    labels <- key$population[ix]
  }
  paste(sort(labels),collapse=";")
}, "")
if (anyNA(hpd$descendants) || any(!nzchar(hpd$descendants)))
  stop("MCC interval extraction produced an empty clade.")
if (anyDuplicated(hpd$descendants))
  stop("Duplicate HPD clades after parsing: ",
       paste(unique(hpd$descendants[duplicated(hpd$descendants)]), collapse=" | "))
hit <- match(node_table$descendants,hpd$descendants)
if (anyNA(hit)) stop("HPD clades differ from this tree; do not transfer intervals across topologies.")
if (any(abs(node_table$age_Ma-hpd$mean_Ma[hit])>1e-5))
  stop("Tree ages differ from the uploaded MCC interval source. Use its pruned all-chromosome tree.")
node_table$lower_Ma <- hpd$lower_Ma[hit]
node_table$upper_Ma <- hpd$upper_Ma[hit]
stopifnot(all(node_table$lower_Ma>=0),all(node_table$upper_Ma<=MAX_AGE),
          all(node_table$lower_Ma<=node_table$upper_Ma))
write.csv(node_table, paste0(OUTPUT_PREFIX, "_node_ages.csv"), row.names=FALSE)
integrifolia <- c("intg_CAswGL", "intg_Nslope", "intg_nGL")
mono <- if (all(integrifolia %in% tr$tip.label)) is.monophyletic(tr, integrifolia) else NA


# Rotate child order toward geography; branch lengths and clades are unchanged.
pop_x <- tapply(d$x[d$in_s170], d$population[d$in_s170], mean)
leaf_order <- function(node) {
  if (node <= n) return(node)
  ch <- children[[as.character(node)]]
  if (REORDER_FOR_GEOGRAPHY) ch <- ch[order(vapply(ch, function(z) mean(pop_x[descendants(z)]), 0.0))]
  unlist(lapply(ch, leaf_order), use.names=FALSE)
}
ordered_tips <- leaf_order(root)
xnode <- rep(NA_real_, n+tr$Nnode)
xnode[ordered_tips] <- seq(-0.95, 0.95, length.out=n)
position <- function(node) {
  if (node <= n) return(xnode[node])
  xnode[node] <<- mean(vapply(children[[as.character(node)]], position, 0.0))
  xnode[node]
}
position(root)
ynode <- TIP_BASELINE + TREE_HEIGHT * ages / MAX_AGE
populations <- sort(tr$tip.label)
species <- sort(unique(unname(pop_species)))
bold <- c("#0072B2", "#D55E00", "#009E73", "#CC0077", "#6F2DA8",
          "#A67600", "#E31A1C", "#00A6B2", "#3B4CC0", "#7A4B24", "#65A30D",
          "#1B365D", "#B34D94", "#005A32", "#FF7F00")
if (length(species)>length(bold)) stop("Extend the species palette.")
species_colors <- setNames(bold[seq_along(species)],species)
cols <- setNames(unname(species_colors[pop_species[populations]]),populations)
d$color <- ifelse(d$in_s170, cols[d$population], "#8D9499")
write.csv(data.frame(population=populations, species=unname(pop_species[populations]), color=unname(cols)), paste0(OUTPUT_PREFIX,"_colors.csv"), row.names=FALSE)
write.csv(d, paste0(OUTPUT_PREFIX,"_sample_audit.csv"), row.names=FALSE)

# The maps package includes coastline data: no separate GIS download is needed.
old_s2 <- sf_use_s2(FALSE)

world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
world <- st_transform(world, 4326)

# Shift land by exactly the same amount as the sample points.
shifted_land <- st_geometry(world) - c(CENTRAL_LONGITUDE, 0)
shifted_land <- st_set_crs(shifted_land, 4326)

# Split polygons at the new map edges before projecting.
shifted_land <- st_wrap_dateline(
  shifted_land,
  options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
  quiet = TRUE
)

world <- st_sf(geometry = shifted_land)
world <- st_make_valid(world)
world <- st_transform(world, crs_map)

# Clip against the existing projected northern-hemisphere boundary.
world <- suppressWarnings(st_intersection(world, boundary))

sf_use_s2(old_s2)

world_geometry <- (st_geometry(world) - c(0, ybottom)) / radius

draw <- function() {
  layout(
  matrix(c(1, 2), ncol = 1),
  heights = c(8.5, 1.5)
)
  par(mar=c(1.2, 0.6, 0.6, 0.15), family="sans", xaxs="i", yaxs="i")
  plot.new(); plot.window(xlim=c(-1.30,1.20), ylim=c(-0.10,TIP_BASELINE+TREE_HEIGHT+0.08), asp=1)
  plot(outline,add=TRUE,col="white",border=NA)
  plot(world_geometry, add=TRUE, col="#E5E5E5", border="#B7B7B7", lwd=0.4)
  plot(outline,add=TRUE,col=NA,border="#4D555B",lwd=0.6)
  ticks <- seq(0,MAX_AGE,by=2)
  for (a in ticks) {
    y <- TIP_BASELINE+TREE_HEIGHT*a/MAX_AGE
    segments(-0.95,y,0.95,y,col="#E3E7E9",lwd=0.6)
    text(-1.065,y,format(a,trim=TRUE),adj=1,cex=0.95,col="#526169")
  }
  text(-1.07,TIP_BASELINE+TREE_HEIGHT+0.055,"Ma",adj=1,cex=1.05,font=2)
  # Per-sample connections, including coincident samples. No locality averaging.
  foreground <- which(d$in_s170)
  for (i in foreground) {
    tip <- match(d$population[i],tr$tip.label)
    segments(xnode[tip],TIP_BASELINE-0.19,d$x[i],d$y[i],col=adjustcolor(d$color[i],alpha.f=CONNECTION_ALPHA),lwd=0.85)
  }
  # Larger hollow rings ensure S380-only records remain visible at shared localities.
  bg <- !d$in_s170
  points(d$x[bg],d$y[bg],pch=1,col="#8D9499",cex=0.78,lwd=0.8)
  points(d$x[foreground],d$y[foreground],pch=21,bg=d$color[foreground],col="#FFFFFF",cex=0.69,lwd=0.35)
  # Age uncertainty is vertical because age is the y coordinate.
  rect(xnode[node_table$node]-0.010,
       TIP_BASELINE+TREE_HEIGHT*node_table$lower_Ma/MAX_AGE,
       xnode[node_table$node]+0.010,
       TIP_BASELINE+TREE_HEIGHT*node_table$upper_Ma/MAX_AGE,
       col=adjustcolor("#527A99",alpha.f=0.35),border=NA)
  for (parent in unique(tr$edge[,1])) {
    ch <- children[[as.character(parent)]]
    segments(min(xnode[ch]),ynode[parent],max(xnode[ch]),ynode[parent],lwd=1,col="#343C40")
    for (child in ch) segments(xnode[child],ynode[child],xnode[child],ynode[parent],lwd=1,
      col=if(child<=n) cols[tr$tip.label[child]] else "#343C40")
  }
  points(xnode[seq_len(n)],rep(TIP_BASELINE,n),pch=16,col=cols[tr$tip.label],cex=0.48)
  # Display names: old population identifier = new printed label
label_changes <- c(
  "octo_Pobeda" = "sp_Pobeda",
  "punc_CHU" = "gran_Chukotka",
  "punc_Yak" = "gran_Yakutia",
  "drum" = "drummondii",
  "cauc" = "caucasica"
)

tip_labels <- gsub("_", " ", tr$tip.label)
i <- match(tr$tip.label, names(label_changes))
tip_labels[!is.na(i)] <- unname(label_changes[i[!is.na(i)]])

text(
  xnode[seq_len(n)], rep(TIP_BASELINE - 0.013, n),
  labels = tip_labels,
  srt = 90, adj = 1, cex = 0.8,
  col = cols[tr$tip.label]
)
  #text(0,-0.065,sprintf("S170: %d samples | Other S380: %d | Display jitter: up to %g km",
   # sum(d$in_s170),sum(bg),JITTER_KM),cex=0.85,col="#444444")
# Full-width bottom legend in six columns.
par(mar = c(0.2, 0.6, 0.2, 0.6), xpd = NA)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

legend(
  "center",
  legend = c(species, "Additional samples"),
  title = "Species",
  text.font = c(rep(3, length(species)), 1),
  title.font = 1,
  pch = c(rep(21, length(species)), 1),
  pt.bg = c(unname(species_colors[species]), NA),
  col = c(rep("white", length(species)), "#8D9499"),
  pt.cex = 1.35,
  bty = "n",
  ncol = 6,
  cex = 1.02,
  x.intersp = 0.65,
  y.intersp = 1.25
)

layout(1)

}
pdf(paste0(OUTPUT_PREFIX,".pdf"),width=16,height=10,useDingbats=FALSE)
draw(); dev.off()
png(paste0(OUTPUT_PREFIX,".png"),width=4800,height=3000,res=300)
draw(); dev.off()
if (interactive()) draw()
message(sprintf("Done: %d samples (%d S170); root %.4f Ma; integrifolia monophyletic: %s",nrow(d),sum(d$in_s170),root_age,mono))
