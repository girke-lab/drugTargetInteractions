## =====================================================================
##  generate_gtopdb_example.R
##  Maintainer-run script that (re)generates the shipped GtoPdb extract
##  under inst/extdata/:
##    gtopdb_interactions_example.tsv
##    gtopdb_hgnc_mapping_example.tsv
##
##  These stand in for GtoPdb's two bulk download files in the man-page
##  examples and as the vignette's offline fallback, so the whole
##  build-then-query path runs without network access. See
##  R/gtoPdbAccess.R's header for why the package builds from the bulk
##  files rather than the (now API-key-gated) REST endpoint, and why
##  shipping an extract is allowed at all - GtoPdb is ODbL with contents
##  under CC-BY-SA 4.0.
##
##  Run manually, NOT as part of any automated build - whenever GtoPdb
##  publishes a release whose shape has changed, or to refresh the
##  extract. Run from the package root:
##    source("inst/scripts/generate_gtopdb_example.R")
##
##  The extract keeps only the columns .gtoPdbFromTsv() reads, which is
##  what keeps it a few KB rather than a few hundred; each file's
##  original "# GtoPdb Version:" first line is preserved verbatim so
##  .gtoPdbVersion() still works on it.
## =====================================================================

outDir <- "inst/extdata"
stopifnot(dir.exists(outDir))

## Genes to keep. FGFR1 and NLRP3 are the two of the project's 8 test
## genes GtoPdb actually covers; pemigatinib's rows bring in its
## FGFR2/FGFR3 targets, which the drug -> target examples assert on.
genes   <- c("FGFR1", "NLRP3")
ligands <- c("pemigatinib")

## Columns .gtoPdbFromTsv() reads, in the order the real file has them.
cols <- c("Target", "Target ID", "Target Species", "Ligand ID", "Ligand",
          "Type", "Action", "Selectivity", "Primary Target",
          "Affinity Units", "Affinity High", "Affinity Median",
          "Affinity Low", "PubMed ID")

base <- "https://www.guidetopharmacology.org/DATA"
ixFile  <- tempfile(fileext = ".tsv")
mapFile <- tempfile(fileext = ".tsv")
utils::download.file(file.path(base, "interactions.tsv"), ixFile, quiet = TRUE)
utils::download.file(file.path(base, "GtP_to_HGNC_mapping.tsv"), mapFile,
                     quiet = TRUE)

readBulk <- function(path) {
    read.delim(path, skip = 1L, sep = "\t", quote = "\"",
               stringsAsFactors = FALSE, check.names = FALSE,
               colClasses = "character")
}
ix  <- readBulk(ixFile)
map <- readBulk(mapFile)
ixLines  <- readLines(ixFile, warn = FALSE, encoding = "UTF-8")
mapLines <- readLines(mapFile, warn = FALSE, encoding = "UTF-8")

## Target IDs for the chosen genes, via the mapping file's target/object
## rows only - its "IUPHAR ID" column also carries ligand IDs in the
## same numeric space (see R/gtoPdbAccess.R's header).
isObject <- grepl("ObjectDisplayForward", map[["GtP URL"]], fixed = TRUE)
tids <- unique(map[["IUPHAR ID"]][isObject & map[["HGNC Symbol"]] %in% genes])

keep <- ix[["Target Species"]] == "Human" &
        nzchar(ix[["Target ID"]]) &
        (ix[["Target ID"]] %in% tids | ix[["Ligand"]] %in% ligands)
stopifnot(sum(keep) > 0L)

quoteField <- function(x) {
    ifelse(is.na(x), "", paste0("\"", gsub("\"", "\"\"", x), "\""))
}
sub <- ix[keep, cols, drop = FALSE]
body <- apply(sub, 1L, function(r) paste(quoteField(r), collapse = "\t"))
writeLines(c(ixLines[1], paste(quoteField(cols), collapse = "\t"), body),
           file.path(outDir, "gtopdb_interactions_example.tsv"),
           useBytes = TRUE)

## Every target/object mapping row for a target ID present in the
## extract, so multi-subunit complexes keep all their symbols.
mapIdx <- which(isObject & map[["IUPHAR ID"]] %in% unique(ix[["Target ID"]][keep]))
stopifnot(length(mapIdx) > 0L)
writeLines(c(mapLines[1], mapLines[2], mapLines[2 + mapIdx]),
           file.path(outDir, "gtopdb_hgnc_mapping_example.tsv"),
           useBytes = TRUE)

message("wrote ", sum(keep), " interaction rows and ", length(mapIdx),
        " mapping rows from GtoPdb ",
        sub(".*Version:\\s*", "", sub(" - published.*", "", ixLines[1])))
