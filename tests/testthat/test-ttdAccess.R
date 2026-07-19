## Network-guarded tests for the TTD flat-file download / local-SQLite
## build / query path. All network-touching tests skip cleanly when
## offline (BiocCheck penalises unguarded network calls). Defined locally
## (not shared from test-apiAccess.R) so this file runs standalone too.

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release - see test-apiAccess.R's copy of this helper
## for the full rationale. Remove selectively later, not all at once.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## Builds (or reuses a cached) TTD SQLite once per test run and shares it
## across every test_that() below, rather than re-downloading/rebuilding
## per test - buildTtdDb() is not free the way a single REST GET is.
.ttdTestDb <- new.env(parent = emptyenv())
.getTtdTestDbPath <- function() {
    if (is.null(.ttdTestDb$path)) .ttdTestDb$path <- buildTtdDb()
    .ttdTestDb$path
}

test_that(".ttdVersion parses the 'Version X.Y.Z' header line, no network", {
    f <- tempfile()
    writeLines(c("TTD flat file", "some legend line",
                "Version 10.1.01 (2024.01.10)", "more legend"), f)
    expect_identical(.ttdVersion(f), "10.1.01")

    f2 <- tempfile()
    writeLines(c("no version line here"), f2)
    expect_true(is.na(.ttdVersion(f2)))
})

test_that(".ttdReadLong/.ttdFieldLookup/.ttdParseTargets parse a synthetic P1-01-style file, no network", {
    f <- tempfile()
    writeLines(c(
        "Legend line 1",
        "Legend line 2",
        "Version 10.1.01 (2024.01.10)",
        "T00001\tGENENAME\tFGFR1",
        "T00001\tUNIPROID\tFGFR1_HUMAN",
        "T00001\tTARGTYPE\tSuccessful",
        "T00001\tDRUGINFO\tD00001\tPemigatinib\tApproved",
        "T00002\tGENENAME\tPTGS1",
        "T00002\tUNIPROID\tPTGS1_HUMAN",
        "T00002\tTARGTYPE\tSuccessful"
    ), f)

    tgt <- .ttdParseTargets(f)
    expect_equal(nrow(tgt), 2L)
    expect_identical(tgt$GeneName[tgt$TargetID == "T00001"], "FGFR1")
    expect_identical(tgt$Uniprot[tgt$TargetID == "T00001"], "FGFR1_HUMAN")

    drugNames <- .ttdParseTargetDrugNames(f)
    expect_identical(drugNames$DrugName[drugNames$DrugID == "D00001"], "Pemigatinib")
})

test_that(".ttdParseDrugSmiles parses a synthetic P1-02-style file, no network", {
    f <- tempfile()
    writeLines(c(
        "Legend line 1", "Legend line 2", "Version 10.1.01 (2024.01.10)",
        "D00001\tDRUGSMIL\tCCN1C2=C3C=C(NC3=NC=C2CN(C1=O)C4=C(C(=CC(=C4F)OC)OC)F)CN5CCOCC5",
        "D00002\tDRUGSMIL\tCC(=O)OC1=CC=CC=C1C(=O)O"
    ), f)
    sm <- .ttdParseDrugSmiles(f)
    expect_equal(nrow(sm), 2L)
    expect_identical(sm$Smiles[sm$DrugID == "D00002"], "CC(=O)OC1=CC=CC=C1C(=O)O")
})

test_that("downloadTTD fetches/caches the 3 TTD flat files", {
    skip_if_offline_dti()
    paths <- downloadTTD()
    expect_identical(names(paths), c("targets", "drugs", "mapping"))
    expect_true(all(vapply(paths, file.exists, logical(1))))
})

test_that("buildTtdDb builds a queryable local SQLite with the expected table", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    expect_true(file.exists(dbPath))
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_true("ttd_interactions" %in% dbListTables(con))
    cols <- dbListFields(con, "ttd_interactions")
    expect_identical(cols, c("TargetID", "GeneName", "Uniprot", "TargetType",
                             "DrugID", "DrugName", "Smiles", "Highest_status", "MOA"))
})

test_that("ttdTargetAnnot matches the validated 8-gene row count, NA-padding undrugged targets", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    df <- ttdTargetAnnot(list(molType = "protein", idType = "symbol", ids = genes), dbPath)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "TargetID", "GeneName", "Uniprot",
                                  "TargetType", "DrugID", "DrugName", "Smiles",
                                  "Highest_status", "MOA"))
    ## 61 real matches (R_Py_code/ttdAccess.R reference) + 3 NA-padded rows
    ## for genes TTD has zero interactions for (ADIPOR1, ADIPOR2, TFEB).
    expect_equal(sum(!is.na(df$TargetID)), 61L)
    expect_equal(nrow(df), 64L)
    naGenes <- df$QueryIDs[is.na(df$TargetID)]
    expect_setequal(naGenes, c("ADIPOR1", "ADIPOR2", "TFEB"))
})

test_that("ttdTargetAnnot resolves pemigatinib's known FGFR targets (drug -> target)", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    df <- ttdTargetAnnot(list(molType = "cmp", idType = "name", ids = "Pemigatinib"), dbPath)
    expect_equal(nrow(df), 3L)
    expect_true(all(df$QueryIDs == "Pemigatinib"))
    expect_setequal(df$GeneName, c("FGFR1", "FGFR2", "FGFR3"))
    expect_true(all(df$Highest_status == "Approved"))
})

test_that("ttdTargetAnnot supports exact-match native IDs (uniprot mnemonic, ttd_target_id, ttd_drug_id)", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    byMnemonic <- ttdTargetAnnot(list(molType = "protein", idType = "uniprot",
                                      ids = "FGFR1_HUMAN"), dbPath)
    expect_true(all(byMnemonic$GeneName == "FGFR1"))

    byTid <- ttdTargetAnnot(list(molType = "protein", idType = "ttd_target_id",
                                 ids = unique(byMnemonic$TargetID)), dbPath)
    expect_equal(nrow(byTid), nrow(byMnemonic))

    byDid <- ttdTargetAnnot(list(molType = "cmp", idType = "ttd_drug_id",
                                 ids = "D0O6UY"), dbPath)
    expect_true(all(byDid$DrugName == "Pemigatinib"))
})

test_that("ttdTargetAnnot surfaces unmatched query IDs as NA rows, incl. all-unmatched", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    allUnmatched <- ttdTargetAnnot(list(molType = "protein", idType = "symbol",
                                        ids = "NOTAREALGENEXYZ"), dbPath)
    expect_equal(nrow(allUnmatched), 1L)
    expect_true(is.na(allUnmatched$TargetID))

    mixed <- ttdTargetAnnot(list(molType = "cmp", idType = "name",
                                 ids = c("Pemigatinib", "NOTAREALDRUGXYZ")), dbPath)
    naRow <- mixed[mixed$QueryIDs == "NOTAREALDRUGXYZ", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$DrugID))
    expect_equal(sum(mixed$QueryIDs == "Pemigatinib"), 3L)
})

test_that("ttdTargetAnnot rejects unsupported idType / malformed queryBy", {
    skip_if_offline_dti()
    dbPath <- .getTtdTestDbPath()
    expect_error(
        ttdTargetAnnot(list(molType = "protein", idType = "ensembl", ids = "x"), dbPath),
        "must be one of")
    expect_error(
        ttdTargetAnnot(list(molType = "cmp", idType = "name", ids = character(0)), dbPath),
        "need to be populated")
})
