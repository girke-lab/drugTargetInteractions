## Tests for the Guide to PHARMACOLOGY (GtoPdb) bulk-download /
## local-SQLite build / query path. The build-and-query tests run offline
## against the small extract of GtoPdb's two bulk files that ships in
## inst/extdata; only the tests that genuinely exercise the download are
## network-guarded, and those skip cleanly when offline (BiocCheck
## penalises unguarded network calls). The guard helper is defined
## locally (not shared from test-apiAccess.R) so this file runs
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

## The small extract of GtoPdb's two bulk files that ships with the
## package. It covers FGFR1 and NLRP3 plus pemigatinib's FGFR rows -
## exactly the subset these tests assert on - so the build-and-query
## path can be exercised in full without touching the network.
.getGtoPdbExamplePaths <- function() {
    list(
        interactions = system.file("extdata", "gtopdb_interactions_example.tsv",
                                   package = "drugTargetInteractions"),
        hgncMapping = system.file("extdata", "gtopdb_hgnc_mapping_example.tsv",
                                  package = "drugTargetInteractions"))
}

## Builds the SQLite once per test run and shares it across every
## test_that() below rather than rebuilding per test.
.gtoPdbTestDb <- new.env(parent = emptyenv())
.getGtoPdbTestDbPath <- function() {
    if (is.null(.gtoPdbTestDb$path))
        .gtoPdbTestDb$path <- buildGtoPdbDb(paths = .getGtoPdbExamplePaths())
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
    dbPath <- .getGtoPdbTestDbPath()
    expect_true(file.exists(dbPath))
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_true("gtp_interactions" %in% dbListTables(con))
    cols <- dbListFields(con, "gtp_interactions")
    expect_identical(cols, c("target_gene", "targetId", "targetName", "species",
                             "primaryTarget", "ligandId", "ligandName", "type",
                             "action", "affinity", "affinityParameter",
                             "selectivity", "pubmedId"))
    ## species filter applied at build time
    allSpecies <- dbGetQuery(con, "SELECT DISTINCT species FROM gtp_interactions")$species
    expect_identical(allSpecies, "Human")
})

test_that("gtoPdbTargetAnnot matches the validated 8-gene coverage, NA-padding uncovered targets", {
    dbPath <- .getGtoPdbTestDbPath()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    df <- gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol", ids = genes), dbPath)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "target_gene", "targetId", "targetName",
                                  "species", "primaryTarget", "ligandId", "ligandName",
                                  "type", "action", "affinity", "affinityParameter",
                                  "selectivity", "pubmedId"))
    ## Live-validated 2026-07-19: only FGFR1 and NLRP3 are covered by
    ## GtoPdb of these 8 (a more targeted/curated resource than TTD/Broad,
    ## not exhaustive) - the other 6 are NA-padded.
    naGenes <- unique(df$QueryIDs[is.na(df$targetId)])
    expect_setequal(naGenes, c("FGF21", "KLB", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2"))
    expect_gt(sum(df$QueryIDs == "FGFR1" & !is.na(df$targetId)), 0L)
    expect_gt(sum(df$QueryIDs == "NLRP3" & !is.na(df$targetId)), 0L)
})

test_that("gtoPdbTargetAnnot resolves pemigatinib's known FGFR targets (drug -> target)", {
    dbPath <- .getGtoPdbTestDbPath()
    df <- gtoPdbTargetAnnot(list(molType = "cmp", idType = "name", ids = "pemigatinib"), dbPath)
    expect_true(all(df$QueryIDs == "pemigatinib"))
    expect_setequal(df$target_gene, c("FGFR1", "FGFR2", "FGFR3"))
    expect_true(all(df$type == "Inhibitor"))
    expect_false(any(grepl("<", df$targetName, fixed = TRUE)))  ## HTML stripped
})

test_that("gtoPdbTargetAnnot supports exact-match native IDs (gtp_target_id, gtp_ligand_id)", {
    dbPath <- .getGtoPdbTestDbPath()
    byTid <- gtoPdbTargetAnnot(list(molType = "protein", idType = "gtp_target_id",
                                    ids = "1808"), dbPath)
    expect_true(all(byTid$target_gene == "FGFR1"))

    byLid <- gtoPdbTargetAnnot(list(molType = "cmp", idType = "gtp_ligand_id",
                                    ids = "9767"), dbPath)
    expect_true(all(byLid$ligandName == "pemigatinib"))
})

test_that("gtoPdbTargetAnnot surfaces unmatched query IDs as NA rows, incl. all-unmatched", {
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
    dbPath <- .getGtoPdbTestDbPath()
    expect_error(
        gtoPdbTargetAnnot(list(molType = "protein", idType = "uniprot", ids = "x"), dbPath),
        "must be one of")
    expect_error(
        gtoPdbTargetAnnot(list(molType = "cmp", idType = "name", ids = character(0)), dbPath),
        "need to be populated")
})

test_that("buildGtoPdbDb(paths=) builds from local files without touching the cache", {
    tmpCache <- tempfile("dtiCache")
    dir.create(tmpCache)
    old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
    on.exit({
        if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR")
        else Sys.setenv(R_USER_CACHE_DIR = old)
        unlink(tmpCache, recursive = TRUE)
    })
    Sys.setenv(R_USER_CACHE_DIR = tmpCache)

    dbPath <- buildGtoPdbDb(paths = .getGtoPdbExamplePaths())
    expect_true(file.exists(dbPath))
    ## Nothing was written into the cache, so a supplied-paths build can
    ## never shadow the real gtopdb_<version>.db.
    expect_false(any(grepl("gtopdb", list.files(tmpCache, recursive = TRUE))))
})

test_that("buildGtoPdbDb rejects a malformed or non-existent 'paths'", {
    expect_error(buildGtoPdbDb(paths = list(interactions = "a")),
                 "must be a named list")
    expect_error(
        buildGtoPdbDb(paths = list(interactions = tempfile(), hgncMapping = tempfile())),
        "do not exist")
})

test_that("the build drops non-human rows and rows carrying no target identity", {
    ex <- .getGtoPdbExamplePaths()
    lines <- readLines(ex$interactions, warn = FALSE, encoding = "UTF-8")
    ## One mouse row and one row with no "Target ID" at all, both of which
    ## must be gone from the built table.
    body <- strsplit(lines[3], "\t", fixed = TRUE)[[1]]
    mouse <- body; mouse[3] <- "\"Mouse\""
    noTarget <- body; noTarget[2] <- "\"\""
    f <- tempfile(fileext = ".tsv")
    writeLines(c(lines, paste(mouse, collapse = "\t"),
                 paste(noTarget, collapse = "\t")), f)

    dbPath <- buildGtoPdbDb(paths = list(interactions = f,
                                         hgncMapping = ex$hgncMapping))
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_identical(
        dbGetQuery(con, "SELECT DISTINCT species FROM gtp_interactions")$species,
        "Human")
    expect_equal(
        dbGetQuery(con, "SELECT COUNT(*) n FROM gtp_interactions WHERE targetId = ''")$n,
        0L)
})

test_that(".gtoPdbAffinity rebuilds a range or a single value, matching the REST payload", {
    ## Bounds that differ become a range; equal or absent bounds collapse
    ## to the single value. Single-precision noise is reported to one
    ## decimal, as GtoPdb's own API did.
    expect_identical(.gtoPdbAffinity("5.0", "5.25", "5.5"), "5.0 - 5.5")
    expect_identical(.gtoPdbAffinity("9.3", "9.3", "9.3"), "9.3")
    expect_identical(.gtoPdbAffinity("", "9.399999618530273", ""), "9.4")
    expect_true(is.na(.gtoPdbAffinity("", "", "")))
})

test_that("downloadGtoPdb reports a missing or rejected API key for its REST route", {
    expect_error(downloadGtoPdb(useServices = TRUE, apiKey = ""),
                 "require a registered API key")
    expect_error(downloadGtoPdb(useServices = TRUE, apiKey = NA_character_),
                 "require a registered API key")
})

test_that("downloadGtoPdb's REST route reports GtoPdb's 401 as an authentication failure", {
    skip_if_offline_dti()
    ## GtoPdb gates every /services/ endpoint behind a key, so an
    ## obviously invalid one is a stable way to exercise the 401 branch.
    expect_error(
        downloadGtoPdb(useServices = TRUE, apiKey = "not-a-real-api-key"),
        "authentication")
})
