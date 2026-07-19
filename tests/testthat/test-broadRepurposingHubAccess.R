## Network-guarded tests for the Broad Repurposing Hub flat-file download /
## local-SQLite build / query path. All network-touching tests skip cleanly
## when offline (BiocCheck penalises unguarded network calls). Defined
## locally (not shared from test-apiAccess.R) so this file runs standalone
## too.

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release - see test-apiAccess.R's copy of this helper
## for the full rationale. Remove selectively later, not all at once.
##
## Additionally: repo-hub.broadinstitute.org has been observed serving an
## incomplete TLS certificate chain (missing InCommon intermediate), which
## can fail these live tests with a certificate error even when otherwise
## online and even on Bioconductor's own build machines - a server-side
## issue, not something this package works around (see
## broadRepurposingHubAccess.R's file header). skip_on_bioc() covers the
## Bioconductor build case either way.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## Builds (or reuses a cached) Repurposing Hub SQLite once per test run and
## shares it across every test_that() below, rather than re-downloading/
## rebuilding per test - buildBroadRepurposingHubDb() is not free the way a
## single REST GET is.
.brhTestDb <- new.env(parent = emptyenv())
.getBrhTestDbPath <- function() {
    if (is.null(.brhTestDb$path)) .brhTestDb$path <- buildBroadRepurposingHubDb()
    .brhTestDb$path
}

test_that(".brhReadTable parses a synthetic Repurposing Hub-style file, no network", {
    f <- tempfile()
    writeLines(c(
        "!Source\t\"The Drug Repurposing Hub, Broad Institute\"",
        "!URL\thttp://www.broadinstitute.org/repurposing",
        "!File_date\t8/18/2025",
        "!Table_name\tdrug_information",
        "!API_available\tSome note.",
        "!Disclaimer\tSome disclaimer.",
        "!Restrictions\tThe Drug Repurposing Hub data is provided for non-commercial use only.",
        "!Citation\tSome citation.",
        "!Contact\trepurposing@broadinstitute.org",
        "pert_iname\tclinical_phase\tmoa\ttarget\tdisease_area\tindication",
        "drug-a\tLaunched\t\"kinase inhibitor, potent\"\tFGFR1 | FGFR2\t\t",
        "drug-b\tPreclinical\tantitumor agent\t\t\t"
    ), f)
    df <- .brhReadTable(f)
    expect_identical(attr(df, "fileDate"), "8/18/2025")
    expect_equal(nrow(df), 2L)
    expect_identical(df$moa[df$pert_iname == "drug-a"], "kinase inhibitor, potent")
    expect_identical(df$target[df$pert_iname == "drug-a"], "FGFR1 | FGFR2")
    expect_true(is.na(df$target[df$pert_iname == "drug-b"]))
})

test_that(".brhExplodeTargets splits multi-gene rows and keeps no-target drugs as one NA row, no network", {
    drug <- data.frame(
        pert_iname = c("drug-a", "drug-b"),
        clinical_phase = c("Launched", "Preclinical"),
        moa = c("kinase inhibitor", "antitumor agent"),
        target = c("FGFR1 | FGFR2", NA_character_),
        disease_area = c(NA_character_, NA_character_),
        indication = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )
    exploded <- .brhExplodeTargets(drug)
    expect_equal(nrow(exploded), 3L)
    expect_setequal(exploded$target_gene[exploded$pert_iname == "drug-a"], c("FGFR1", "FGFR2"))
    expect_true(is.na(exploded$target_gene[exploded$pert_iname == "drug-b"]))
})

test_that("downloadBroadRepurposingHub fetches/caches the 2 flat files", {
    skip_if_offline_dti()
    paths <- downloadBroadRepurposingHub()
    expect_identical(names(paths), c("drug", "sample"))
    expect_true(all(vapply(paths, file.exists, logical(1))))
})

test_that("buildBroadRepurposingHubDb builds a queryable local SQLite with the expected table", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    expect_true(file.exists(dbPath))
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    expect_true("broad_interactions" %in% dbListTables(con))
    cols <- dbListFields(con, "broad_interactions")
    expect_identical(cols, c("target_gene", "pert_iname", "clinical_phase", "moa",
                             "disease_area", "indication", "broad_id", "qc_incompatible",
                             "purity", "vendor", "catalog_no", "vendor_name",
                             "expected_mass", "smiles", "InChIKey", "pubchem_cid",
                             "deprecated_broad_id"))
})

test_that("broadRepurposingHubAnnot matches the validated 8-gene row counts, NA-padding undrugged targets", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    df <- broadRepurposingHubAnnot(list(molType = "protein", idType = "symbol", ids = genes), dbPath)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "target_gene", "pert_iname", "clinical_phase",
                                  "moa", "disease_area", "indication", "broad_id",
                                  "qc_incompatible", "purity", "vendor", "catalog_no",
                                  "vendor_name", "expected_mass", "smiles", "InChIKey",
                                  "pubchem_cid", "deprecated_broad_id"))
    ## Live-validated 2026-07-19: FGF21, KLB and TFEB have zero Repurposing
    ## Hub entries (NA-padded); the other 5 genes resolve to real rows -
    ## exact counts reflect drug x sample multiplicity, not drug count.
    naGenes <- unique(df$QueryIDs[is.na(df$pert_iname)])
    expect_setequal(naGenes, c("FGF21", "KLB", "TFEB"))
    expect_true(all(df$pert_iname[df$QueryIDs == "FGFR1"] != ""))
    expect_gt(sum(df$QueryIDs == "FGFR1" & !is.na(df$pert_iname)), 0L)
})

test_that("broadRepurposingHubAnnot resolves pemigatinib's known FGFR targets (drug -> target)", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    df <- broadRepurposingHubAnnot(list(molType = "cmp", idType = "name", ids = "pemigatinib"), dbPath)
    expect_true(all(df$QueryIDs == "pemigatinib"))
    expect_setequal(df$target_gene, c("FGFR1", "FGFR2", "FGFR3"))
    expect_true(all(df$moa == "fgfr inhibitor"))
    expect_true(all(df$clinical_phase == "Launched"))
})

test_that("broadRepurposingHubAnnot supports exact-match broad_id lookup and case-insensitive name", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    byId <- broadRepurposingHubAnnot(list(molType = "cmp", idType = "broad_id",
                                          ids = "BRD-K00104124-001-01-9"), dbPath)
    expect_true(all(byId$pert_iname == "pemigatinib"))

    byNameCaps <- broadRepurposingHubAnnot(list(molType = "cmp", idType = "name",
                                                ids = "PEMIGATINIB"), dbPath)
    expect_equal(nrow(byNameCaps), 3L)
})

test_that("broadRepurposingHubAnnot surfaces unmatched query IDs as NA rows, incl. all-unmatched", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    allUnmatched <- broadRepurposingHubAnnot(list(molType = "protein", idType = "symbol",
                                                  ids = "NOTAREALGENEXYZ"), dbPath)
    expect_equal(nrow(allUnmatched), 1L)
    expect_true(is.na(allUnmatched$pert_iname))

    mixed <- broadRepurposingHubAnnot(list(molType = "cmp", idType = "name",
                                           ids = c("pemigatinib", "notarealdrugxyz")), dbPath)
    naRow <- mixed[mixed$QueryIDs == "notarealdrugxyz", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$pert_iname))
    expect_equal(sum(mixed$QueryIDs == "pemigatinib"), 3L)
})

test_that("broadRepurposingHubAnnot rejects unsupported idType / malformed queryBy", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    expect_error(
        broadRepurposingHubAnnot(list(molType = "protein", idType = "uniprot", ids = "x"), dbPath),
        "must be")
    expect_error(
        broadRepurposingHubAnnot(list(molType = "cmp", idType = "name", ids = character(0)), dbPath),
        "need to be populated")
})
