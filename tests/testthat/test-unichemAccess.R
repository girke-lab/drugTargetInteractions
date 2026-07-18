## Tests for the UniChem local-SQLite compound cross-referencing layer.
## Deliberately does NOT call buildUnichemDb()/downloadUnichemTables()
## live: the real reference.tsv.gz is ~1.5GB / ~166M rows and takes real
## time to download and process (live-validated manually, see
## PROGRESS.md), far too heavy to run routinely in a test suite. Instead:
##  - .unichemLoadReference()/.unichemFilterAndIndex() are tested against
##    a tiny synthetic reference.tsv.gz fixture built in-test.
##  - getUnichemMapping() has no network dependency at all (it only ever
##    touches a local SQLite path), so it is tested end-to-end against a
##    small hand-built fixture database - no skip_if_offline_dti() needed
##    anywhere in this file.

test_that(".unichemLatestCachedDb finds the most recently added unichem_*.db regardless of its date suffix", {
    bfc <- .getCache()
    f1 <- tempfile(fileext = ".db"); file.create(f1)
    f2 <- tempfile(fileext = ".db"); file.create(f2)
    rid1 <- names(bfcadd(bfc, "unichem_TEST_OLD.db", f1, action = "copy"))
    Sys.sleep(1.1)  ## ensure a distinct create_time from rid1
    rid2 <- names(bfcadd(bfc, "unichem_TEST_NEW.db", f2, action = "copy"))
    on.exit(bfcremove(bfc, c(rid1, rid2)), add = TRUE)

    latest <- .unichemLatestCachedDb()
    expect_identical(basename(latest), basename(bfcrpath(bfc, rids = rid2)))
})

.buildFixtureDb <- function() {
    dbPath <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    xref <- data.frame(
        uci        = c(1L, 1L, 1L,   2L, 2L,   3L),
        src_id     = c(1L, 22L, 50L, 1L, 2L,   22L),
        compound_id = c("CHEMBL25", "2244", "AIN",
                        "CHEMBL1421", "DB01254",
                        "999999"),
        stringsAsFactors = FALSE
    )
    dbWriteTable(con, "unichem_xref", xref, overwrite = TRUE)
    sources <- data.frame(
        src_id   = c(1L, 2L, 22L, 50L),
        name     = c("chembl", "drugbank", "pubchem", "CCDC"),  ## CCDC deliberately uppercase
        name_long = c("ChEMBL", "DrugBank", "PubChem Compounds", "CSD"),
        description = c("d1", "d2", "d3", "d4"),
        stringsAsFactors = FALSE
    )
    dbWriteTable(con, "unichem_sources", sources, overwrite = TRUE)
    dbPath
}

test_that("getUnichemMapping translates via a shared UCI, source names case-insensitive", {
    dbPath <- .buildFixtureDb()
    ## aspirin: chembl -> pubchem
    r1 <- getUnichemMapping("CHEMBL25", from = "chembl", to = "pubchem", dbPath)
    expect_identical(r1$From, "CHEMBL25")
    expect_identical(r1$To, "2244")

    ## chembl -> CCDC, whose stored name is uppercase - query lowercase
    r2 <- getUnichemMapping("CHEMBL25", from = "chembl", to = "ccdc", dbPath)
    expect_identical(r2$To, "AIN")
    ## and uppercase query for the same, matching the stored casing exactly
    r3 <- getUnichemMapping("CHEMBL25", from = "CHEMBL", to = "CCDC", dbPath)
    expect_identical(r3$To, "AIN")
})

test_that("getUnichemMapping accepts numeric source IDs interchangeably with names", {
    dbPath <- .buildFixtureDb()
    byName <- getUnichemMapping("CHEMBL1421", from = "chembl", to = "drugbank", dbPath)
    byId   <- getUnichemMapping("CHEMBL1421", from = 1, to = 2, dbPath)
    expect_identical(byName, byId)
    expect_identical(byId$To, "DB01254")
})

test_that("getUnichemMapping returns no row for an unmapped or unknown compound", {
    dbPath <- .buildFixtureDb()
    r <- getUnichemMapping("NOTAREALID", from = "chembl", to = "pubchem", dbPath)
    expect_equal(nrow(r), 0L)
    expect_identical(names(r), c("From", "To"))
})

test_that("getUnichemMapping errors clearly on an unrecognised source name", {
    dbPath <- .buildFixtureDb()
    expect_error(
        getUnichemMapping("CHEMBL25", from = "chembl", to = "notarealsource", dbPath),
        "Unrecognised UniChem source")
})

test_that("getUnichemMapping handles multiple input ids in one call", {
    dbPath <- .buildFixtureDb()
    ## "2244" (uci=1) has a chembl cross-ref; "999999" (uci=3) does not
    r <- getUnichemMapping(c("2244", "999999"), from = "pubchem", to = "chembl", dbPath)
    expect_equal(nrow(r), 1L)
    expect_identical(r$From, "2244")
    expect_identical(r$To, "CHEMBL25")
})

test_that(".unichemLoadReference chunk-loads a synthetic reference.tsv.gz correctly", {
    f <- tempfile(fileext = ".tsv.gz")
    con <- gzcon(file(f, "wb"))
    writeLines(c(
        "UCI\tSRC_ID\tSRC_COMPOUND_ID\tASSIGMENT",
        "1\t1\tCHEMBL25\t1",
        "1\t22\t2244\t1",
        "2\t41\tSLM:000627379\t0"  ## non-numeric compound_id, e.g. SwissLipids
    ), con)
    close(con)

    dbPath <- tempfile(fileext = ".db")
    dbcon <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(dbcon))
    n <- .unichemLoadReference(f, dbcon, chunkSize = 2L)  ## force multi-chunk
    expect_equal(n, 3L)
    got <- dbGetQuery(dbcon, "SELECT * FROM ref_raw ORDER BY uci, src_id")
    expect_equal(nrow(got), 3L)
    expect_identical(got$compound_id[got$src_id == 41], "SLM:000627379")
})

test_that(".unichemFilterAndIndex keeps only multi-source UCIs and drops the staging table", {
    dbPath <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    raw <- data.frame(
        uci = c(1L, 1L, 2L, 3L, 3L, 3L),
        src_id = c(1L, 22L, 1L, 1L, 2L, 22L),
        compound_id = c("A", "B", "C", "D", "E", "F"),
        stringsAsFactors = FALSE
    )
    dbWriteTable(con, "ref_raw", raw, overwrite = TRUE)
    .unichemFilterAndIndex(con, minSources = 2L, anchorSources = NULL)

    expect_false("ref_raw" %in% dbListTables(con))
    expect_true("unichem_xref" %in% dbListTables(con))
    kept <- dbGetQuery(con, "SELECT DISTINCT uci FROM unichem_xref ORDER BY uci")$uci
    expect_identical(kept, c(1L, 3L))  ## uci=2 (single source) dropped
    expect_equal(dbGetQuery(con, "SELECT COUNT(*) AS n FROM unichem_xref")$n, 5L)
})

test_that(".unichemFilterAndIndex's anchorSources further restricts to drug-relevant UCIs", {
    ## Mirrors the real-world case found live 2026-07-17: a UCI can have
    ## >=2 sources (satisfying minSources) while none of them is
    ## drug-relevant (e.g. pubchem+surechembl only) - anchorSources
    ## must exclude that UCI even though minSources alone would keep it.
    dbPath <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    raw <- data.frame(
        uci = c(1L, 1L,   2L, 2L,   3L, 3L, 3L),
        src_id = c(1L, 22L,  15L, 22L,  1L, 22L, 15L),
        compound_id = c("CHEMBL1", "PC1",  "SC2", "PC2",  "CHEMBL3", "PC3", "SC3"),
        stringsAsFactors = FALSE
    )
    dbWriteTable(con, "ref_raw", raw, overwrite = TRUE)
    sources <- data.frame(
        src_id = c(1L, 15L, 22L),
        name = c("chembl", "surechembl", "pubchem"),
        name_long = c("ChEMBL", "SureChEMBL", "PubChem Compounds"),
        description = c("d1", "d2", "d3"),
        stringsAsFactors = FALSE
    )
    dbWriteTable(con, "unichem_sources", sources, overwrite = TRUE)

    .unichemFilterAndIndex(con, minSources = 2L, anchorSources = "chembl")

    kept <- dbGetQuery(con, "SELECT DISTINCT uci FROM unichem_xref ORDER BY uci")$uci
    ## uci=1 and uci=3 both have a chembl row; uci=2 (surechembl+pubchem
    ## only, no chembl) satisfies minSources but not the anchor and must
    ## be dropped.
    expect_identical(kept, c(1L, 3L))
})

test_that(".unichemResolveSourceId is case-insensitive and rejects unknown names", {
    dbPath <- .buildFixtureDb()
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_equal(.unichemResolveSourceId("chembl", con), 1L)
    expect_equal(.unichemResolveSourceId("CHEMBL", con), 1L)
    expect_equal(.unichemResolveSourceId("ccdc", con), 50L)
    expect_equal(.unichemResolveSourceId("22", con), 22L)  ## numeric passthrough, no lookup
    expect_error(.unichemResolveSourceId("notarealsource", con), "Unrecognised UniChem source")
})
