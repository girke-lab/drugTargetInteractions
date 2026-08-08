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

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release - see test-apiAccess.R's copy of this helper
## for the full rationale. Remove selectively later, not all at once.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
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

.brhTestDb <- new.env(parent = emptyenv())
.getBrhTestDbPath <- function() {
    if (is.null(.brhTestDb$path)) .brhTestDb$path <- buildBroadRepurposingHubDb()
    .brhTestDb$path
}

.gtoPdbTestDb <- new.env(parent = emptyenv())
.getGtoPdbTestDbPath <- function() {
    if (is.null(.gtoPdbTestDb$path)) .gtoPdbTestDb$path <- buildGtoPdbDb()
    .gtoPdbTestDb$path
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
    expect_identical(.dtiMetaQueryBy("dgidb", FALSE, "aspirin"),
                     list(molType = "cmp", idType = "name", ids = "aspirin"))
    expect_identical(.dtiMetaQueryBy("opentargets", FALSE, "CHEMBL25"),
                     list(molType = "cmp", idType = "name", ids = "CHEMBL25"))
    expect_identical(.dtiMetaQueryBy("broad", TRUE, "FGFR1"),
                     list(molType = "protein", idType = "symbol", ids = "FGFR1"))
    expect_identical(.dtiMetaQueryBy("broad", FALSE, "pemigatinib"),
                     list(molType = "cmp", idType = "name", ids = "pemigatinib"))
    expect_identical(.dtiMetaQueryBy("gtopdb", TRUE, "FGFR1"),
                     list(molType = "protein", idType = "symbol", ids = "FGFR1"))
    expect_identical(.dtiMetaQueryBy("gtopdb", FALSE, "pemigatinib"),
                     list(molType = "cmp", idType = "name", ids = "pemigatinib"))
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
                            sources = c("chembl", "dgidb", "opentargets"))
    expect_true(setequal(names(res), c("chembl", "dgidb", "opentargets")))
    expect_identical(nrow(res$dgidb), nrow(getDgidbDrugs("FGFR1")))
    expect_identical(nrow(res$opentargets), nrow(getOpenTargetsDrugs("FGFR1")))
    resolved <- attr(res, "resolved")
    expect_identical(unname(resolved$chembl), "P11362")
    expect_identical(unname(resolved$dgidb), "FGFR1")
})

test_that("queryDrugTargets rejects 'pubchem' as a source (bioassay data, not an annotation source)", {
    expect_error(
        queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                         sources = "pubchem"),
        "sources")
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

test_that("queryDrugTargets dispatches to Broad Repurposing Hub given a local db path, matching broadRepurposingHubAnnot directly", {
    skip_if_offline_dti()
    dbPath <- .getBrhTestDbPath()
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                            sources = "broad", brhDbPath = dbPath)
    direct <- broadRepurposingHubAnnot(list(molType = "protein", idType = "symbol", ids = "FGFR1"), dbPath)
    expect_identical(res$broad, direct)
})

test_that("queryDrugTargets dispatches to GtoPdb given a local db path, matching gtoPdbTargetAnnot directly", {
    skip_if_offline_dti()
    dbPath <- .getGtoPdbTestDbPath()
    res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
                            sources = "gtopdb", gtoPdbDbPath = dbPath)
    direct <- gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol", ids = "FGFR1"), dbPath)
    expect_identical(res$gtopdb, direct)
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

## --- combineDrugTargets() --------------------------------------------------

## Synthetic queryDrugTargets()-shaped fixture, network-free: real column
## names for every source, tiny data. Mimics starting from an Ensembl ID
## that got resolved differently per source (uniprot for chembl, symbol
## for the rest) - exercises the "resolved" attribute's original-query-id
## backfill, not just a same-string passthrough.
.combineFixture <- function() {
    res <- list(
        chembl = data.frame(
            QueryIDs = "P11362", chembl_id = "CHEMBL1421", Drug_Name = "DASATINIB",
            MOA = "Tyrosine-protein kinase ABL inhibitor", Action_Type = "INHIBITOR",
            Max_Phase = 4, First_Approval = 2006, ChEMBL_TID = "CHEMBL1862",
            UniProt_ID = "P11362", Desc = "d", Organism = "Homo sapiens",
            Mesh_Indication = NA, stringsAsFactors = FALSE),
        dgidb = data.frame(
            QueryIDs = "FGFR1", gene_name = "FGFR1", drug_name = "SUNITINIB",
            drug_concept_id = "chembl:CHEMBL535", drug_approved = TRUE,
            interaction_types = "inhibitor", directionality = "INHIBITORY",
            interaction_score = 1.2, evidence_score = 3.4, sources = "ChEMBL",
            db = "DGIdb", stringsAsFactors = FALSE),
        opentargets = data.frame(
            QueryIDs = "FGFR1", ensembl_id = "ENSG00000077782", approved_symbol = "FGFR1",
            drug_id = "CHEMBL1201733", drug_name = "PAZOPANIB", drug_type = "Small molecule",
            max_clinical_stage = 4, mechanism_of_action = "FGFR inhibitor",
            action_type = "INHIBITOR", disease_id = NA, disease_name = NA,
            stringsAsFactors = FALSE),
        ttd = data.frame(
            QueryIDs = "FGFR1", TargetID = "T47101", GeneName = "FGFR1",
            Uniprot = "FGFR1_HUMAN", TargetType = "Successful", DrugID = "D00ABO",
            DrugName = "KW-2449", Smiles = "x", Highest_status = "Phase 1", MOA = "Inhibitor",
            stringsAsFactors = FALSE),
        broad = data.frame(
            QueryIDs = "FGFR1", target_gene = "FGFR1", pert_iname = "pemigatinib",
            clinical_phase = "Launched", moa = "fgfr inhibitor", disease_area = NA,
            indication = NA, smiles = "x", InChIKey = "x", pubchem_cid = "1",
            structure_ambiguous = FALSE, stringsAsFactors = FALSE),
        gtopdb = data.frame(
            QueryIDs = "FGFR1", target_gene = "FGFR1", targetId = 1808L,
            targetName = "fibroblast growth factor receptor 1", species = "Human",
            primaryTarget = TRUE, ligandId = 9767L, ligandName = "pemigatinib",
            type = "Inhibitor", action = "Inhibition", affinity = "7.0",
            affinityParameter = "pIC50", selectivity = NA, refIds = "34085",
            stringsAsFactors = FALSE)
    )
    attr(res, "resolved") <- list(
        chembl = c(ENSG00000077782 = "P11362"),
        dgidb = c(ENSG00000077782 = "FGFR1"),
        opentargets = c(ENSG00000077782 = "FGFR1"),
        ttd = c(ENSG00000077782 = "FGFR1"),
        broad = c(ENSG00000077782 = "FGFR1"),
        gtopdb = c(ENSG00000077782 = "FGFR1")
    )
    res
}

test_that("combineDrugTargets maps every source's columns correctly and backfills the original query_id", {
    combined <- combineDrugTargets(.combineFixture())
    expect_equal(nrow(combined), 6L)
    expect_identical(unique(combined$query_id), "ENSG00000077782")  ## not each source's own resolved id
    expect_identical(combined$gene_symbol[combined$source == "ChEMBL"], NA_character_)
    expect_identical(combined$gene_symbol[combined$source == "DGIdb"], "FGFR1")
    expect_identical(combined$drug_name[combined$source == "ChEMBL"], "DASATINIB")
    expect_identical(combined$drug_name[combined$source == "TTD"], "KW-2449")
    expect_identical(combined$drug_name[combined$source == "Broad Repurposing Hub"], "pemigatinib")
    expect_identical(combined$drug_name[combined$source == "GtoPdb"], "pemigatinib")
    expect_identical(combined$action[combined$source == "DGIdb"], "inhibitor")
    expect_identical(combined$action[combined$source == "Broad Repurposing Hub"], "fgfr inhibitor")
    expect_identical(combined$action[combined$source == "GtoPdb"], "Inhibition")
    expect_setequal(combined$source, c("ChEMBL", "DGIdb", "OpenTargets", "TTD",
                                       "Broad Repurposing Hub", "GtoPdb"))
})

test_that("combineDrugTargets rejects 'pubchem' as a results name (bioassay data, not combinable)", {
    fixture <- .combineFixture()
    fixture$pubchem <- fixture$dgidb
    expect_error(combineDrugTargets(fixture), "does not recognise source")
})

test_that("combineDrugTargets falls back to each source's own QueryIDs when there is no 'resolved' attribute", {
    fixture <- .combineFixture()
    attr(fixture, "resolved") <- NULL
    combined <- combineDrugTargets(fixture)
    expect_identical(combined$query_id[combined$source == "ChEMBL"], "P11362")
    expect_identical(combined$query_id[combined$source == "DGIdb"], "FGFR1")
})

test_that("combineDrugTargets respects a custom columns subset", {
    combined <- combineDrugTargets(.combineFixture(), columns = c("drug_name", "source"))
    expect_identical(names(combined), c("drug_name", "source"))
})

test_that("combineDrugTargets(resolveGeneSymbol=TRUE) fills ChEMBL's gene_symbol without touching other sources", {
    skip_if_offline_dti()
    combined <- combineDrugTargets(.combineFixture(), resolveGeneSymbol = TRUE)
    expect_identical(combined$gene_symbol[combined$source == "ChEMBL"], "FGFR1")
    expect_identical(combined$gene_symbol[combined$source == "DGIdb"], "FGFR1")
})

test_that("combineDrugTargets returns an empty, correctly-columned data.frame for an empty results list", {
    empty <- combineDrugTargets(list())
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), c("query_id", "gene_symbol", "drug_name", "action", "source"))
})

test_that("combineDrugTargets errors clearly on an unrecognised source name", {
    fixture <- .combineFixture()
    names(fixture)[1] <- "notarealsource"
    expect_error(combineDrugTargets(fixture), "does not recognise source")
})

test_that("combineDrugTargets errors clearly when an expected column is missing from a source's data", {
    fixture <- .combineFixture()
    fixture$chembl$Drug_Name <- NULL
    expect_error(combineDrugTargets(fixture), "expected column 'Drug_Name' not found")
})

test_that("combineDrugTargets tolerates a source whose IDs resolved but that matched zero rows", {
    ## Regression guard: queryDrugTargets() can hand back a non-NULL,
    ## 0-row data.frame for a source (e.g. a transient upstream-API gap
    ## after resolution succeeded) - combineDrugTargets() must skip it,
    ## not crash on $<-.data.frame's "replacement has 1 row, data has 0"
    ## when adding the internal .source column to a 0-row frame.
    fixture <- .combineFixture()
    fixture$chembl <- fixture$chembl[0, ]
    combined <- combineDrugTargets(fixture)
    expect_false("ChEMBL" %in% combined$source)
    expect_equal(nrow(combined), 5L)
    expect_setequal(combined$source, c("DGIdb", "OpenTargets", "TTD",
                                       "Broad Repurposing Hub", "GtoPdb"))
})

test_that("combineDrugTargets returns an empty, correctly-columned data.frame when every source matched zero rows", {
    fixture <- .combineFixture()
    fixture <- lapply(fixture, function(df) df[0, ])
    attr(fixture, "resolved") <- attr(.combineFixture(), "resolved")
    empty <- combineDrugTargets(fixture)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), c("query_id", "gene_symbol", "drug_name", "action", "source"))
})
