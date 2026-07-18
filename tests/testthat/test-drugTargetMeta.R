## Tests for the cross-source dispatcher (queryDrugTargets() and its
## gene/compound ID-resolution helpers). Network-guarded throughout,
## following the shared skip_if_offline_dti() pattern - defined locally
## so this file runs standalone.
##
## Compound-side tests that need a UniChem SQLite (structured-ID-to-
## structured-ID / structured-ID-to-name resolution) look for an already
## -built one via BiocFileCache and skip if none is cached - building one
## from scratch takes on the order of an hour (see test-unichemAccess.R),
## far too heavy to trigger from a test suite. TTD is cheap enough to
## build on demand (mirrors test-ttdAccess.R's memoized helper).

skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

.getCachedUnichemDb <- function() {
    bfc <- .getCache()
    hits <- bfcquery(bfc, "unichem_", field = "rname")
    if (nrow(hits) == 0L) return(NA_character_)
    ## most recently added
    hits <- hits[order(hits$create_time, decreasing = TRUE), ]
    tryCatch(.getCacheFile(hits$rname[1]), error = function(e) NA_character_)
}

skip_if_no_unichem_db <- function() {
    testthat::skip_on_cran()
    if (is.na(.getCachedUnichemDb()))
        testthat::skip("No cached UniChem SQLite - build one with buildUnichemDb() to run this test")
}

.ttdTestDb <- new.env(parent = emptyenv())
.getTtdTestDbPath <- function() {
    if (is.null(.ttdTestDb$path)) .ttdTestDb$path <- buildTtdDb()
    .ttdTestDb$path
}

## --- .resolveGeneIds() ---------------------------------------------------

test_that(".resolveGeneIds passthrough needs no network", {
    r <- .resolveGeneIds(c("FGFR1", "KLB"), idType = "symbol", to = "symbol")
    expect_identical(unname(r), c("FGFR1", "KLB"))
    expect_identical(names(r), c("FGFR1", "KLB"))
})

test_that(".resolveGeneIds resolves symbol -> uniprot (1-hop), NA for unmapped", {
    skip_if_offline_dti()
    r <- .resolveGeneIds(c("FGFR1", "NOTAREALGENEXYZ"), idType = "symbol", to = "uniprot")
    expect_identical(unname(r["FGFR1"]), "P11362")
    expect_true(is.na(r["NOTAREALGENEXYZ"]))
})

test_that(".resolveGeneIds routes symbol <-> ensembl through uniprot (2-hop)", {
    skip_if_offline_dti()
    toEnsembl <- .resolveGeneIds("FGFR1", idType = "symbol", to = "ensembl")
    expect_identical(unname(toEnsembl), "ENSG00000077782.24")
    toSymbol <- .resolveGeneIds("ENSG00000077782", idType = "ensembl", to = "symbol")
    expect_identical(unname(toSymbol), "FGFR1")
})

test_that(".resolveGeneIds rejects an unrecognised idType/to", {
    expect_error(.resolveGeneIds("FGFR1", idType = "notatype", to = "symbol"))
    expect_error(.resolveGeneIds("FGFR1", idType = "symbol", to = "notatype"))
})

## --- .resolveCompoundIds() ------------------------------------------------

test_that(".resolveCompoundIds passthrough needs no network", {
    r <- .resolveCompoundIds("CHEMBL25", idType = "chembl_id", to = "chembl_id")
    expect_identical(unname(r), "CHEMBL25")
})

test_that(".resolveCompoundIds structured -> name via getChemblMolecule, needs no UniChem", {
    skip_if_offline_dti()
    r <- .resolveCompoundIds("CHEMBL25", idType = "chembl_id", to = "name")
    expect_identical(unname(r), "ASPIRIN")
})

test_that(".resolveCompoundIds name -> pubchem_id needs no UniChem (direct PubChem resolver)", {
    skip_if_offline_dti()
    r <- .resolveCompoundIds("aspirin", idType = "name", to = "pubchem_id")
    expect_identical(unname(r), "2244")
})

test_that(".resolveCompoundIds errors clearly when a structured-to-structured hop needs UniChem and none is supplied", {
    expect_error(
        .resolveCompoundIds("CHEMBL25", idType = "chembl_id", to = "pubchem_id"),
        "requires unichemDbPath")
})

test_that(".resolveCompoundIds structured -> structured via UniChem, and 2-hop structured -> name", {
    skip_if_offline_dti()
    skip_if_no_unichem_db()
    dbPath <- .getCachedUnichemDb()
    toPubchem <- .resolveCompoundIds("CHEMBL25", idType = "chembl_id", to = "pubchem_id",
                                     unichemDbPath = dbPath)
    expect_identical(unname(toPubchem), "2244")
    toDrugbank <- .resolveCompoundIds("CHEMBL25", idType = "chembl_id", to = "drugbank_id",
                                      unichemDbPath = dbPath)
    expect_identical(unname(toDrugbank), "DB00945")
    ## 2-hop: pubchem_id -> chembl_id -> name
    toName <- .resolveCompoundIds("2244", idType = "pubchem_id", to = "name",
                                  unichemDbPath = dbPath)
    expect_identical(unname(toName), "ASPIRIN")
})

## --- .dtiMetaQueryBy() (pure, no network) ---------------------------------

test_that(".dtiMetaQueryBy builds the correct native queryBy per source/direction", {
    expect_identical(.dtiMetaQueryBy("chembl", TRUE, "P11362"),
                     list(molType = "protein", idType = "Uniprot", ids = "P11362"))
    expect_identical(.dtiMetaQueryBy("ttd", TRUE, "FGFR1"),
                     list(molType = "protein", idType = "symbol", ids = "FGFR1"))
    expect_identical(.dtiMetaQueryBy("pubchem", FALSE, "aspirin"),
                     list(molType = "cmp", idType = "name", ids = "aspirin"))
    expect_identical(.dtiMetaQueryBy("opentargets", FALSE, "CHEMBL25"),
                     list(molType = "cmp", idType = "name", ids = "CHEMBL25"))
})

## --- queryDrugTargets() ----------------------------------------------------

test_that("queryDrugTargets rejects malformed queryBy without any network calls", {
    expect_error(
        queryDrugTargets(list(molType = "gene", idType = "symbol", ids = character(0))),
        "need to be populated")
    expect_error(
        queryDrugTargets(list(molType = "nonsense", idType = "symbol", ids = "FGFR1")),
        "molType")
    expect_error(
        queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                         sources = "notarealsource"),
        "sources")
})

test_that("queryDrugTargets: gene direction dispatches to multiple sources and matches their standalone counts", {
    skip_if_offline_dti()
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                            sources = c("chembl", "pubchem", "dgidb", "opentargets"))
    expect_true(setequal(names(res), c("chembl", "pubchem", "dgidb", "opentargets")))
    expect_identical(nrow(res$pubchem), nrow(getPubchemDrugs("FGFR1")))
    expect_identical(nrow(res$dgidb), nrow(getDgidbDrugs("FGFR1")))
    expect_identical(nrow(res$opentargets), nrow(getOpenTargetsDrugs("FGFR1")))
    resolved <- attr(res, "resolved")
    expect_identical(unname(resolved$chembl), "P11362")
    expect_identical(unname(resolved$pubchem), "FGFR1")
})

test_that("queryDrugTargets resolves a non-native starting ID type (Ensembl) via 2-hop routing", {
    skip_if_offline_dti()
    res <- queryDrugTargets(list(molType = "gene", idType = "ensembl", ids = "ENSG00000077782"),
                            sources = "opentargets")
    expect_equal(nrow(res$opentargets), 94L)
    expect_identical(unname(attr(res, "resolved")$opentargets), "FGFR1")
})

test_that("queryDrugTargets: compound direction with a UniChem-resolved ChEMBL ID matches getChemblDrugTarget", {
    skip_if_offline_dti()
    skip_if_no_unichem_db()
    dbPath <- .getCachedUnichemDb()
    res <- queryDrugTargets(list(molType = "cmp", idType = "drugbank_id", ids = "DB00945"),
                            sources = "chembl", unichemDbPath = dbPath)
    direct <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id", ids = "CHEMBL25"))
    expect_identical(res$chembl, direct)
})

test_that("queryDrugTargets dispatches to TTD given a local TTD db path, matching ttdTargetAnnot directly", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                            sources = "ttd", ttdDbPath = dbPath)
    direct <- ttdTargetAnnot(list(molType = "protein", idType = "symbol", ids = "FGFR1"), dbPath)
    expect_identical(res$ttd, direct)
})

test_that("queryDrugTargets returns an empty list (not an error) when a required local db path is missing", {
    ## symbol->symbol is a zero-network passthrough, so this needs no
    ## internet access - the point is purely the missing-ttdDbPath path.
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                            sources = "ttd")
    expect_identical(res, structure(list(), resolved = list(ttd = c(FGFR1 = "FGFR1"))))
})

test_that("queryDrugTargets returns an empty list when nothing resolves for the requested source", {
    skip_if_offline_dti()
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "NOTAREALGENEXYZ"),
                            sources = "chembl")
    expect_length(res, 0L)
})
