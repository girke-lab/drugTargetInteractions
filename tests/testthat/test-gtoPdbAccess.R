## Network-guarded tests for the Guide to PHARMACOLOGY (GtoPdb) bulk-REST
## download / local-SQLite build / query path. All network-touching tests
## skip cleanly when offline (BiocCheck penalises unguarded network calls).
## Defined locally (not shared from test-apiAccess.R) so this file runs
## standalone too.

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release - see test-apiAccess.R's copy of this helper
## for the full rationale. Remove selectively later, not all at once.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## Builds (or reuses a cached) GtoPdb SQLite once per test run and shares
## it across every test_that() below, rather than re-downloading/
## rebuilding per test - buildGtoPdbDb() is not free the way a single
## REST GET is (the bulk interactions pull is ~12MB).
.gtoPdbTestDb <- new.env(parent = emptyenv())
.getGtoPdbTestDbPath <- function() {
    if (is.null(.gtoPdbTestDb$path)) .gtoPdbTestDb$path <- buildGtoPdbDb()
    .gtoPdbTestDb$path
}

test_that(".gtoPdbStripHtml removes markup while keeping inner text, no network", {
    expect_identical(.gtoPdbStripHtml("5-HT<sub>1A</sub> receptor"), "5-HT1A receptor")
    expect_identical(.gtoPdbStripHtml("plain text"), "plain text")
    expect_identical(.gtoPdbStripHtml("<i>x</i> <small>y</small>"), "x y")
})

test_that(".gtoPdbVersion parses the '# GtoPdb Version: X.Y' header line, no network", {
    f <- tempfile()
    writeLines(c("# GtoPdb Version: 2026.2 - published: 2026-06-15",
                "\"HGNC Symbol\"\t\"HGNC ID\""), f)
    expect_identical(.gtoPdbVersion(f), "2026.2")

    f2 <- tempfile()
    writeLines("no version line here", f2)
    expect_true(is.na(.gtoPdbVersion(f2)))
})

test_that(".gtoPdbReadHgncMapping filters to target rows only, dropping ligand-ID collisions, no network", {
    f <- tempfile()
    writeLines(c(
        "# GtoPdb Version: 2026.2 - published: 2026-06-15",
        "\"HGNC Symbol\"\t\"HGNC ID\"\t\"IUPHAR Name\"\t\"IUPHAR ID\"\t\"GtP URL\"",
        "\"ASIC1\"\t\"100\"\t\"ASIC1\"\t\"684\"\thttps://x/GRAC/ObjectDisplayForward?objectId=684",
        "\"ADM2\"\t\"28898\"\t\"adrenomedullin 2\"\t\"684\"\thttps://x/GRAC/LigandDisplayForward?ligandId=684",
        "\"ATP5MC1\"\t\"841\"\t\"F-type ATPase C subunit\"\t\"803\"\thttps://x/GRAC/ObjectDisplayForward?objectId=803",
        "\"ATP5MC2\"\t\"842\"\t\"F-type ATPase C subunit\"\t\"803\"\thttps://x/GRAC/ObjectDisplayForward?objectId=803"
    ), f)
    map <- .gtoPdbReadHgncMapping(f)
    expect_equal(nrow(map), 3L)
    expect_identical(map$target_gene[map$targetId == 684], "ASIC1")  ## not "ADM2", the ligand-ID collision
    expect_setequal(map$target_gene[map$targetId == 803], c("ATP5MC1", "ATP5MC2"))
})

test_that("downloadGtoPdb fetches/caches the bulk interactions + HGNC mapping files", {
    skip_if_offline_dti()
    paths <- downloadGtoPdb()
    expect_identical(names(paths), c("interactions", "hgncMapping"))
    expect_true(all(vapply(paths, file.exists, logical(1))))
})

test_that("buildGtoPdbDb builds a queryable local SQLite with the expected table", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    expect_true(file.exists(dbPath))
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_true("gtp_interactions" %in% dbListTables(con))
    cols <- dbListFields(con, "gtp_interactions")
    expect_identical(cols, c("target_gene", "targetId", "targetName", "species",
                             "primaryTarget", "ligandId", "ligandName", "type",
                             "action", "affinity", "affinityParameter",
                             "selectivity", "refIds"))
    ## species filter applied at build time
    allSpecies <- dbGetQuery(con, "SELECT DISTINCT species FROM gtp_interactions")$species
    expect_identical(allSpecies, "Human")
})

test_that("gtoPdbTargetAnnot matches the validated 8-gene coverage, NA-padding uncovered targets", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    df <- gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol", ids = genes), dbPath)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "target_gene", "targetId", "targetName",
                                  "species", "primaryTarget", "ligandId", "ligandName",
                                  "type", "action", "affinity", "affinityParameter",
                                  "selectivity", "refIds"))
    ## Live-validated 2026-07-19: only FGFR1 and NLRP3 are covered by
    ## GtoPdb of these 8 (a more targeted/curated resource than TTD/Broad,
    ## not exhaustive) - the other 6 are NA-padded.
    naGenes <- unique(df$QueryIDs[is.na(df$targetId)])
    expect_setequal(naGenes, c("FGF21", "KLB", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2"))
    expect_gt(sum(df$QueryIDs == "FGFR1" & !is.na(df$targetId)), 0L)
    expect_gt(sum(df$QueryIDs == "NLRP3" & !is.na(df$targetId)), 0L)
})

test_that("gtoPdbTargetAnnot resolves pemigatinib's known FGFR targets (drug -> target)", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    df <- gtoPdbTargetAnnot(list(molType = "cmp", idType = "name", ids = "pemigatinib"), dbPath)
    expect_true(all(df$QueryIDs == "pemigatinib"))
    expect_setequal(df$target_gene, c("FGFR1", "FGFR2", "FGFR3"))
    expect_true(all(df$type == "Inhibitor"))
    expect_false(any(grepl("<", df$targetName, fixed = TRUE)))  ## HTML stripped
})

test_that("gtoPdbTargetAnnot supports exact-match native IDs (gtp_target_id, gtp_ligand_id)", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    byTid <- gtoPdbTargetAnnot(list(molType = "protein", idType = "gtp_target_id",
                                    ids = "1808"), dbPath)
    expect_true(all(byTid$target_gene == "FGFR1"))

    byLid <- gtoPdbTargetAnnot(list(molType = "cmp", idType = "gtp_ligand_id",
                                    ids = "9767"), dbPath)
    expect_true(all(byLid$ligandName == "pemigatinib"))
})

test_that("gtoPdbTargetAnnot surfaces unmatched query IDs as NA rows, incl. all-unmatched", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    allUnmatched <- gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol",
                                           ids = "NOTAREALGENEXYZ"), dbPath)
    expect_equal(nrow(allUnmatched), 1L)
    expect_true(is.na(allUnmatched$targetId))

    mixed <- gtoPdbTargetAnnot(list(molType = "cmp", idType = "name",
                                    ids = c("pemigatinib", "notarealdrugxyz")), dbPath)
    naRow <- mixed[mixed$QueryIDs == "notarealdrugxyz", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$ligandId))
    expect_equal(sum(mixed$QueryIDs == "pemigatinib"), 3L)
})

test_that("gtoPdbTargetAnnot rejects unsupported idType / malformed queryBy", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    expect_error(
        gtoPdbTargetAnnot(list(molType = "protein", idType = "uniprot", ids = "x"), dbPath),
        "must be one of")
    expect_error(
        gtoPdbTargetAnnot(list(molType = "cmp", idType = "name", ids = character(0)), dbPath),
        "need to be populated")
})
