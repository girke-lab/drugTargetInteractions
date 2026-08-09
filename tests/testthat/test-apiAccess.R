## Network-guarded tests for the live ChEMBL REST accessors. All tests
## skip cleanly when offline so the Bioconductor build machines never
## fail on a flaky endpoint (BiocCheck penalises unguarded network calls).

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release: skip every live-network test on Bioconductor's
## own build machines for now, rather than risk a flaky external API
## (rate limiting, a brief outage) failing the official build. Once
## we've seen clean Bioconductor build reports, remove skip_on_bioc()
## (not this whole function) from individual tests selectively, source
## by source, as each proves reliable there - not all at once.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## For getChemblDrugTarget()'s unichemDbPath = ... path: looks for an
## already-built UniChem SQLite via BiocFileCache and skips if none is
## cached - building one from scratch takes on the order of an hour (see
## test-unichemAccess.R), far too heavy to trigger from this suite.
.getCachedUnichemDb <- function() {
    bfc <- .getCache()
    hits <- bfcquery(bfc, "unichem_", field = "rname")
    if (nrow(hits) == 0L) return(NA_character_)
    hits <- hits[order(hits$create_time, decreasing = TRUE), ]
    tryCatch(.getCacheFile(hits$rname[1]), error = function(e) NA_character_)
}

skip_if_no_unichem_db <- function() {
    testthat::skip_on_cran()
    if (is.na(.getCachedUnichemDb()))
        testthat::skip("No cached UniChem SQLite - build one with buildUnichemDb() to run this test")
}

test_that("getChemblMolecule returns tidy rows with expected columns", {
    skip_if_offline_dti()
    df <- getChemblMolecule(c("CHEMBL25", "CHEMBL1201585"))
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 2L)
    expect_true(all(c("chembl_id", "pref_name", "canonical_smiles",
                      "mw_freebase", "qed_weighted") %in% names(df)))
    expect_identical(df$pref_name[df$chembl_id == "CHEMBL25"], "ASPIRIN")
})

test_that("getChemblTarget resolves a gene symbol to target IDs", {
    skip_if_offline_dti()
    tg <- getChemblTarget("FGFR1", organism = "Homo sapiens")
    expect_s3_class(tg, "data.frame")
    expect_true(nrow(tg) >= 1L)
    expect_true(any(tg$target_type == "SINGLE PROTEIN"))
})

test_that("getChemblBioactivities honours the row cap and type filter", {
    skip_if_offline_dti()
    act <- getChemblBioactivities("CHEMBL1862", standardType = "IC50",
                                  maxRows = 25L)
    expect_s3_class(act, "data.frame")
    expect_lte(nrow(act), 25L)
    if (nrow(act) > 0L) expect_true(all(act$standard_type == "IC50"))
})

test_that("getChemblDrugTarget target->drug returns FGFR1 inhibitors", {
    skip_if_offline_dti()
    df <- getChemblDrugTarget(list(molType = "protein", idType = "Uniprot",
                                   ids = "P11362"))
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "chembl_id", "Drug_Name",
                                  "PubChem_CID", "MOA", "Action_Type",
                                  "Max_Phase", "First_Approval",
                                  "ChEMBL_TID", "UniProt_ID", "Desc",
                                  "Organism", "Mesh_Indication"))
    expect_true(all(is.na(df$PubChem_CID)))  # no unichemDbPath passed
    expect_true(nrow(df) >= 1L)
    expect_true(all(df$QueryIDs == "P11362"))
    expect_true("CHEMBL1201733" %in% df$chembl_id)  # pazopanib
})

test_that("getChemblDrugTarget(unichemDbPath = ...) resolves PubChem_CID via UniChem", {
    skip_if_offline_dti()
    skip_if_no_unichem_db()
    df <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                   ids = "CHEMBL25"),  # aspirin
                              unichemDbPath = .getCachedUnichemDb())
    expect_identical(df$PubChem_CID[df$chembl_id == "CHEMBL25"][1], "2244")
})

test_that("getChemblDrugTarget drug->target returns dasatinib's targets", {
    skip_if_offline_dti()
    df <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                   ids = "CHEMBL1421"))
    expect_s3_class(df, "data.frame")
    expect_true(nrow(df) >= 1L)
    expect_true(all(df$QueryIDs == "CHEMBL1421"))
    expect_true("P00519" %in% df$UniProt_ID)  # ABL1
})

test_that("getChemblDrugTarget surfaces unmatched query IDs as NA rows", {
    skip_if_offline_dti()
    df <- getChemblDrugTarget(list(molType = "protein", idType = "Uniprot",
                                   ids = c("P11362", "NOTAREALACCESSION")))
    expect_true("NOTAREALACCESSION" %in% df$QueryIDs)
    naRow <- df[df$QueryIDs == "NOTAREALACCESSION", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$chembl_id))
})

test_that("getChemblDrugTarget rejects unsupported idType / malformed queryBy", {
    expect_error(
        getChemblDrugTarget(list(molType = "cmp", idType = "PubChem_ID",
                                 ids = "2244")),
        "currently supports only")
    expect_error(
        getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                 ids = character(0))),
        "need to be populated")
})

test_that(".dtiChunk splits IDs into batches of at most `size`, no network", {
    x <- LETTERS[1:7]
    ch <- .dtiChunk(x, 3L)
    expect_equal(unname(lengths(ch)), c(3L, 3L, 1L))
    expect_identical(unlist(ch, use.names = FALSE), x)
    expect_identical(.dtiChunk(character(0), 3L), list())
})

test_that("getChemblDrugTarget batches multiple target IDs identically to one-at-a-time calls", {
    skip_if_offline_dti()
    ## Note: real ChEMBL fan-out is large here (3 accessions -> ~100+ drug
    ## rows), so this deliberately uses the default chunkSize (one request
    ## per lookup stage) rather than a small one -- a tiny chunkSize against
    ## this much fan-out would mean dozens of live requests per stage.
    ## Multi-chunk pagination itself is covered network-cheaply below by
    ## the getChemblMolecule chunkSize=1 case (1 record per ID, no fan-out).
    ids <- c("P11362", "P00519", "P00533")  # FGFR1, ABL1, EGFR
    batched <- getChemblDrugTarget(list(molType = "protein", idType = "Uniprot",
                                        ids = ids))
    expect_true(setequal(unique(batched$QueryIDs), ids))
    single <- getChemblDrugTarget(list(molType = "protein", idType = "Uniprot",
                                       ids = "P11362"))
    key <- function(d) sort(paste(d$chembl_id, d$ChEMBL_TID))
    expect_identical(key(batched[batched$QueryIDs == "P11362", ]), key(single))
})

test_that("getChemblDrugTarget batches multiple compound IDs (drug->target)", {
    skip_if_offline_dti()
    ids <- c("CHEMBL1421", "CHEMBL25")  # dasatinib, aspirin
    df <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                   ids = ids))
    expect_true(setequal(unique(df$QueryIDs), ids))
    expect_true("P00519" %in% df$UniProt_ID[df$QueryIDs == "CHEMBL1421"])
})

test_that("getChemblBioassay target->drug returns FGFR1 IC50 measurements", {
    skip_if_offline_dti()
    df <- getChemblBioassay(list(molType = "protein", idType = "Uniprot",
                                 ids = "P11362"), standardType = "IC50")
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("QueryIDs", "chembl_id", "Drug_Name", "ChEMBL_TID",
                                  "UniProt_ID", "Organism", "Desc", "assay_chembl_id",
                                  "assay_description", "standard_type",
                                  "standard_relation", "standard_value",
                                  "standard_units", "pchembl_value"))
    expect_true(nrow(df) > 0L)
    expect_true(all(df$QueryIDs == "P11362"))
    expect_true(all(df$standard_type == "IC50"))
})

test_that("getChemblBioassay drug->target returns dasatinib's measurements", {
    skip_if_offline_dti()
    df <- getChemblBioassay(list(molType = "cmp", idType = "chembl_id",
                                 ids = "CHEMBL1421"))
    expect_true(nrow(df) > 0L)
    expect_true(all(df$QueryIDs == "CHEMBL1421"))
    expect_identical(unique(df$chembl_id), "CHEMBL1421")
})

test_that("getChemblBioassay surfaces unmatched query IDs as NA rows", {
    skip_if_offline_dti()
    df <- getChemblBioassay(list(molType = "protein", idType = "Uniprot",
                                 ids = c("P11362", "NOTAREALACCESSION")),
                            standardType = "IC50")
    expect_true("NOTAREALACCESSION" %in% df$QueryIDs)
    naRow <- df[df$QueryIDs == "NOTAREALACCESSION", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$chembl_id))
})

test_that("getChemblBioassay rejects unsupported idType / malformed queryBy", {
    expect_error(
        getChemblBioassay(list(molType = "cmp", idType = "PubChem_ID", ids = "2244")),
        "currently supports only")
    expect_error(
        getChemblBioassay(list(molType = "cmp", idType = "chembl_id", ids = character(0))),
        "need to be populated")
})

test_that("getChemblBioassay(fields = 'all') adds activity.*-prefixed columns", {
    skip_if_offline_dti()
    ## One live call only: fields = "core" is by construction a fixed
    ## subset of fields = "all" (.dtiSelectFields()), and the core column
    ## list is already documented with zero network calls via
    ## listBioassayFields() - no need to re-fetch a "core" baseline live
    ## (that path is already covered by the dedicated core test above).
    all <- getChemblBioassay(list(molType = "cmp", idType = "chembl_id", ids = "CHEMBL1421"),
                             fields = "all")
    coreCols <- c("QueryIDs", "chembl_id", "Drug_Name", "ChEMBL_TID", "UniProt_ID",
                 "Organism", "Desc", "assay_chembl_id", "assay_description",
                 "standard_type", "standard_relation", "standard_value",
                 "standard_units", "pchembl_value")
    expect_identical(coreCols, listBioassayFields("chembl")[seq_along(coreCols)])
    expect_true(all(coreCols %in% names(all)))
    expect_true(ncol(all) > length(coreCols))
    expect_true("activity.bao_label" %in% names(all))

    ## Custom-vector selection is client-side logic on the already-
    ## fetched wide data - test it directly, no second live call needed.
    sub <- .dtiSelectFields(all, c("Drug_Name", "activity.assay_type"), coreCols)
    expect_identical(names(sub), c("QueryIDs", "Drug_Name", "activity.assay_type"))
})

test_that("listBioassayFields returns the documented static column list, no network", {
    fc <- listBioassayFields("chembl")
    expect_true("standard_type" %in% fc)
    fp <- listBioassayFields("pubchem")
    expect_true("gene_symbol" %in% fp)
    expect_error(listBioassayFields("dgidb"))
})

test_that("getChemblMolecule batches across chunks, preserves order/duplicates, NA-fills unresolved IDs", {
    skip_if_offline_dti()
    ## molecule.json is 1 record per ID (no fan-out), so chunkSize=1 here
    ## is only 3 cheap requests -- this is what actually exercises the
    ## multi-chunk pagination path end-to-end against the live API.
    df <- getChemblMolecule(c("CHEMBL25", "NOTAREALID", "CHEMBL1421", "CHEMBL25"),
                            chunkSize = 1L)
    expect_equal(nrow(df), 4L)
    expect_identical(df$chembl_id, c("CHEMBL25", "NOTAREALID", "CHEMBL1421", "CHEMBL25"))
    expect_true(is.na(df$pref_name[df$chembl_id == "NOTAREALID"][1]))
    expect_identical(df$pref_name[df$chembl_id == "CHEMBL25"], c("ASPIRIN", "ASPIRIN"))
})

test_that("getPubchemDrugs returns FGFR1 bioactivities incl. AZD4547", {
    skip_if_offline_dti()
    df <- getPubchemDrugs("FGFR1")
    expect_s3_class(df, "data.frame")
    expect_true(all(c("gene_symbol", "geneid", "target_accession", "cid",
                      "drug_name", "canonical_smiles", "activity_name",
                      "activity_value_uM", "assay_name", "db") %in% names(df)))
    expect_true(nrow(df) >= 1L)
    expect_true(all(df$gene_symbol == "FGFR1"))
    expect_identical(unique(df$geneid), "2260")
    expect_true(any(grepl("4547", df$drug_name)))  # AZD4547
    expect_true(all(df$db == "PubChem"))
})

test_that("getPubchemTargets resolves aspirin's human targets, filtered by taxid", {
    skip_if_offline_dti()
    df <- getPubchemTargets("aspirin")
    expect_s3_class(df, "data.frame")
    expect_true(all(c("cid", "drug_name", "gene_symbol", "geneid", "taxid",
                      "target_accession", "activity_name", "activity_value_uM",
                      "assay_name", "db") %in% names(df)))
    expect_true(nrow(df) >= 1L)
    expect_true(all(df$cid == "2244"))
    expect_true(all(df$taxid == 9606L))
    expect_true(all(c("PTGS1", "PTGS2") %in% df$gene_symbol))  # COX1/COX2
})

test_that("getPubchemTargets(taxid = NULL) does not filter by species", {
    skip_if_offline_dti()
    all_sp <- getPubchemTargets("aspirin", taxid = NULL)
    human  <- getPubchemTargets("aspirin", taxid = 9606L)
    expect_true(nrow(all_sp) >= nrow(human))
})

test_that("getPubchemDrugTarget target->drug matches getPubchemDrugs and tags QueryIDs", {
    skip_if_offline_dti()
    df <- getPubchemDrugTarget(list(molType = "gene", idType = "symbol",
                                    ids = "FGFR1"))
    expect_s3_class(df, "data.frame")
    expect_identical(names(df)[1], "QueryIDs")
    expect_true(all(df$QueryIDs == "FGFR1"))
    plain <- getPubchemDrugs("FGFR1")
    expect_equal(nrow(df), nrow(plain))
})

test_that("getPubchemDrugTarget drug->target matches getPubchemTargets and tags QueryIDs", {
    skip_if_offline_dti()
    df <- getPubchemDrugTarget(list(molType = "cmp", idType = "name",
                                    ids = "aspirin"))
    expect_true(all(df$QueryIDs == "aspirin"))
    plain <- getPubchemTargets("aspirin")
    expect_equal(nrow(df), nrow(plain))
    expect_true(all(c("PTGS1", "PTGS2") %in% df$gene_symbol))
})

test_that("getPubchemDrugTarget surfaces unmatched query IDs as NA rows (both directions)", {
    skip_if_offline_dti()
    genes <- getPubchemDrugTarget(list(molType = "gene", idType = "symbol",
                                       ids = c("KLB", "NOTAREALGENEXYZ")))
    expect_true("NOTAREALGENEXYZ" %in% genes$QueryIDs)
    naRow <- genes[genes$QueryIDs == "NOTAREALGENEXYZ", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$cid))

    drugs <- getPubchemDrugTarget(list(molType = "cmp", idType = "name",
                                       ids = c("aspirin", "NOTAREALDRUGXYZ")))
    expect_true("NOTAREALDRUGXYZ" %in% drugs$QueryIDs)
    naRow2 <- drugs[drugs$QueryIDs == "NOTAREALDRUGXYZ", ]
    expect_equal(nrow(naRow2), 1L)
    expect_true(is.na(naRow2$cid))
})

test_that("getPubchemDrugTarget rejects unsupported idType / malformed queryBy", {
    expect_error(
        getPubchemDrugTarget(list(molType = "cmp", idType = "PubChem_ID",
                                  ids = "2244")),
        "currently supports only")
    expect_error(
        getPubchemDrugTarget(list(molType = "cmp", idType = "name",
                                  ids = character(0))),
        "need to be populated")
})

test_that(".dtiPubchemFilterRows keeps only Active + numeric + potency rows, no network", {
    cols <- list("Activity Outcome", "Activity Name", "Activity Value [uM]",
                "CID", "Target Accession", "Target GeneID", "Assay Name")
    rows <- list(
        list(Cell = list("Active", "IC50", "0.5", "1", "P1", "1", "A1")),
        list(Cell = list("Inactive", "IC50", "0.5", "2", "P1", "1", "A1")),
        list(Cell = list("Active", "Solubility", "0.5", "3", "P1", "1", "A1")),
        list(Cell = list("Active", "Ki", "not-a-number", "4", "P1", "1", "A1")))
    df <- .dtiPubchemFilterRows(cols, rows)
    expect_equal(nrow(df), 1L)
    expect_identical(df$CID, "1")
})

test_that("getDgidbDrugs matches the validated 8-gene row count exactly", {
    skip_if_offline_dti()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    df <- getDgidbDrugs(genes)
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("gene_name", "drug_name", "drug_concept_id",
                                  "drug_aliases", "drug_approved", "interaction_types",
                                  "directionality", "interaction_score",
                                  "evidence_score", "sources", "db"))
    expect_equal(nrow(df), 227L)  # matches R_Py_code/dgidb_fetch.py reference
    expect_true(all(df$db == "DGIdb"))
})

test_that("getDgidbTargets resolves imatinib/aspirin's known targets", {
    skip_if_offline_dti()
    df <- getDgidbTargets(c("imatinib", "aspirin"))
    expect_s3_class(df, "data.frame")
    expect_true(nrow(df) >= 1L)
    expect_true(setequal(unique(df$drug_name), c("IMATINIB", "ASPIRIN")))
    expect_true("ABL1" %in% df$gene_name[df$drug_name == "IMATINIB"])
})

test_that("DGIdb normalizes matched names to canonical casing regardless of query case", {
    skip_if_offline_dti()
    genes <- getDgidbDrugs("fgfr1")
    expect_true(all(genes$gene_name == "FGFR1"))
    drugs <- getDgidbTargets("Aspirin")
    expect_true(all(drugs$drug_name == "ASPIRIN"))
})

test_that("getDgidbDrugTarget target->drug matches getDgidbDrugs and tags QueryIDs case-insensitively", {
    skip_if_offline_dti()
    df <- getDgidbDrugTarget(list(molType = "gene", idType = "symbol",
                                  ids = "fgfr1"))
    expect_s3_class(df, "data.frame")
    expect_identical(names(df)[1], "QueryIDs")
    expect_true(all(df$QueryIDs == "fgfr1"))  # echoes original query casing
    expect_true(all(df$gene_name == "FGFR1"))  # canonical DGIdb casing
    plain <- getDgidbDrugs("fgfr1")
    expect_equal(nrow(df), nrow(plain))
})

test_that("getDgidbDrugTarget drug->target matches getDgidbTargets and tags QueryIDs case-insensitively", {
    skip_if_offline_dti()
    df <- getDgidbDrugTarget(list(molType = "cmp", idType = "name",
                                  ids = "Aspirin"))
    expect_true(all(df$QueryIDs == "Aspirin"))
    expect_true(all(df$drug_name == "ASPIRIN"))
    plain <- getDgidbTargets("Aspirin")
    expect_equal(nrow(df), nrow(plain))
})

test_that("getDgidbDrugTarget surfaces unmatched query IDs as NA rows (both directions)", {
    skip_if_offline_dti()
    genes <- getDgidbDrugTarget(list(molType = "gene", idType = "symbol",
                                     ids = c("KLB", "NOTAREALGENEXYZ")))
    expect_true("NOTAREALGENEXYZ" %in% genes$QueryIDs)
    naRow <- genes[genes$QueryIDs == "NOTAREALGENEXYZ", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$drug_name))

    drugs <- getDgidbDrugTarget(list(molType = "cmp", idType = "name",
                                     ids = c("aspirin", "NOTAREALDRUGXYZ")))
    expect_true("NOTAREALDRUGXYZ" %in% drugs$QueryIDs)
    naRow2 <- drugs[drugs$QueryIDs == "NOTAREALDRUGXYZ", ]
    expect_equal(nrow(naRow2), 1L)
    expect_true(is.na(naRow2$gene_name))
})

test_that("getDgidbDrugTarget rejects unsupported idType / malformed queryBy", {
    expect_error(
        getDgidbDrugTarget(list(molType = "gene", idType = "GeneID",
                                ids = "2260")),
        "currently supports only")
    expect_error(
        getDgidbDrugTarget(list(molType = "cmp", idType = "name",
                                ids = character(0))),
        "need to be populated")
})

test_that(".dtiParseDgidbInteractions parses/collapses raw nodes, no network", {
    nodes <- list(list(
        drug = list(name = "ASPIRIN", conceptId = "chembl:CHEMBL25", approved = TRUE,
                    drugAliases = list(list(alias = "CHEMBL:CHEMBL25"),
                                      list(alias = "DRUGBANK:DB00945"))),
        gene = list(name = "PTGS1"),
        interactionScore = 1.5,
        evidenceScore = 3.2,
        interactionTypes = list(list(type = "inhibitor", directionality = "INHIBITORY")),
        sources = list(list(sourceDbName = "ChEMBL"), list(sourceDbName = "DrugBank"))
    ))
    df <- .dtiParseDgidbInteractions(nodes)
    expect_equal(nrow(df), 1L)
    expect_identical(df$gene_name, "PTGS1")
    expect_identical(df$drug_approved, TRUE)
    expect_identical(df$interaction_types, "inhibitor")
    expect_identical(df$sources, "ChEMBL; DrugBank")
    expect_identical(df$drug_aliases, "CHEMBL:CHEMBL25; DRUGBANK:DB00945")
    expect_identical(.dtiParseDgidbInteractions(NULL), .dtiEmptyDgidb())
})

test_that(".dtiParseDgidbInteractions leaves drug_aliases empty (not erroring) when a drug has none", {
    nodes <- list(list(
        drug = list(name = "OBSCURE COMPOUND", conceptId = "rxcui:1", approved = FALSE),
        gene = list(name = "GENEX"),
        interactionScore = NULL, evidenceScore = NULL,
        interactionTypes = list(), sources = list()
    ))
    df <- .dtiParseDgidbInteractions(nodes)
    expect_identical(df$drug_aliases, "")
})

test_that("getOpenTargetsIds resolves gene symbols to Ensembl IDs", {
    skip_if_offline_dti()
    ids <- getOpenTargetsIds(c("FGFR1", "NOT_A_GENE"))
    expect_identical(names(ids), c("FGFR1", "NOT_A_GENE"))
    expect_identical(unname(ids["FGFR1"]), "ENSG00000077782")
    expect_true(is.na(ids["NOT_A_GENE"]))
})

test_that("getOpenTargetsDrugIds resolves drug names to ChEMBL IDs, CHEMBL passes through", {
    skip_if_offline_dti()
    ids <- getOpenTargetsDrugIds(c("aspirin", "NOT_A_DRUG"))
    expect_identical(unname(ids["aspirin"]), "CHEMBL25")
    expect_true(is.na(ids["NOT_A_DRUG"]))
})

test_that("getOpenTargetsDrugs matches the validated FGFR1 row count exactly", {
    skip_if_offline_dti()
    df <- getOpenTargetsDrugs("FGFR1")
    expect_s3_class(df, "data.frame")
    expect_identical(names(df), c("ensembl_id", "approved_symbol", "drug_id",
                                  "drug_name", "drug_type", "max_clinical_stage",
                                  "mechanism_of_action", "action_type",
                                  "disease_id", "disease_name"))
    expect_equal(nrow(df), 94L)  # matches R_Py_code/apiAccess.R reference
    expect_true(all(df$approved_symbol == "FGFR1"))
})

test_that("getOpenTargetsTargets resolves aspirin's COX1/COX2 targets", {
    skip_if_offline_dti()
    df <- getOpenTargetsTargets("aspirin")
    expect_s3_class(df, "data.frame")
    expect_true(setequal(df$approved_symbol, c("PTGS1", "PTGS2")))
    expect_true(all(df$chembl_id == "CHEMBL25"))
})

test_that("getOpenTargetsDrugTarget target->drug matches getOpenTargetsDrugs and tags QueryIDs", {
    skip_if_offline_dti()
    df <- getOpenTargetsDrugTarget(list(molType = "gene", idType = "symbol",
                                        ids = "FGFR1"))
    expect_s3_class(df, "data.frame")
    expect_identical(names(df)[1], "QueryIDs")
    expect_true(all(df$QueryIDs == "FGFR1"))
    plain <- getOpenTargetsDrugs("FGFR1")
    expect_equal(nrow(df), nrow(plain))
})

test_that("getOpenTargetsDrugTarget drug->target matches getOpenTargetsTargets and tags QueryIDs", {
    skip_if_offline_dti()
    df <- getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
                                        ids = "aspirin"))
    expect_true(all(df$QueryIDs == "aspirin"))
    expect_true(setequal(df$approved_symbol, c("PTGS1", "PTGS2")))
    plain <- getOpenTargetsTargets("aspirin")
    expect_equal(nrow(df), nrow(plain))
})

test_that("getOpenTargetsDrugTarget resolves native IDs (ENSG/CHEMBL) without re-resolving", {
    skip_if_offline_dti()
    df <- getOpenTargetsDrugTarget(list(molType = "gene", idType = "symbol",
                                        ids = "ENSG00000077782"))
    expect_true(all(df$QueryIDs == "ENSG00000077782"))
    expect_true(all(df$ensembl_id == "ENSG00000077782"))
})

test_that("getOpenTargetsDrugTarget surfaces unmatched query IDs as NA rows (both directions)", {
    skip_if_offline_dti()
    genes <- getOpenTargetsDrugTarget(list(molType = "gene", idType = "symbol",
                                           ids = c("FGFR1", "NOTAREALGENEXYZ")))
    expect_true("NOTAREALGENEXYZ" %in% genes$QueryIDs)
    naRow <- genes[genes$QueryIDs == "NOTAREALGENEXYZ", ]
    expect_equal(nrow(naRow), 1L)
    expect_true(is.na(naRow$drug_name))

    drugs <- getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
                                           ids = c("aspirin", "NOTAREALDRUGXYZ")))
    expect_true("NOTAREALDRUGXYZ" %in% drugs$QueryIDs)
    naRow2 <- drugs[drugs$QueryIDs == "NOTAREALDRUGXYZ", ]
    expect_equal(nrow(naRow2), 1L)
    expect_true(is.na(naRow2$approved_symbol))
})

test_that("getOpenTargetsDrugTarget rejects unsupported idType / malformed queryBy", {
    expect_error(
        getOpenTargetsDrugTarget(list(molType = "gene", idType = "ensembl",
                                      ids = "ENSG00000077782")),
        "currently supports only")
    expect_error(
        getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
                                      ids = character(0))),
        "need to be populated")
})

test_that(".dtiParseTargetDrugs/.dtiParseDrugTargets expand rows per `expand`, no network", {
    tgt <- list(
        approvedSymbol = "PTGS1",
        drugAndClinicalCandidates = list(rows = list(list(
            maxClinicalStage = 4,
            drug = list(id = "CHEMBL25", name = "ASPIRIN", drugType = "Small molecule",
                       mechanismsOfAction = list(rows = list(
                           list(mechanismOfAction = "Cyclooxygenase inhibitor",
                               actionType = "INHIBITOR")))),
            diseases = list(list(disease = list(id = "EFO_1", name = "Pain")),
                            list(disease = list(id = "EFO_2", name = "Fever")))
        )))
    )
    byMech <- .dtiParseTargetDrugs(tgt, "ENSG00000095303", "mechanism")
    expect_equal(nrow(byMech), 1L)
    expect_identical(byMech$disease_name, "Pain; Fever")
    byDisease <- .dtiParseTargetDrugs(tgt, "ENSG00000095303", "disease")
    expect_equal(nrow(byDisease), 2L)
    expect_identical(.dtiParseTargetDrugs(NULL, "x", "drug"), .dtiEmptyDrugs())

    drg <- list(name = "ASPIRIN", drugType = "Small molecule", maximumClinicalStage = 4,
               mechanismsOfAction = list(rows = list(list(
                   mechanismOfAction = "Cyclooxygenase inhibitor", actionType = "INHIBITOR",
                   targetName = "Prostaglandin G/H synthase",
                   targets = list(list(id = "ENSG00000095303", approvedSymbol = "PTGS1"),
                                 list(id = "ENSG00000073756", approvedSymbol = "PTGS2"))))))
    byTarget <- .dtiParseDrugTargets(drg, "CHEMBL25", "target")
    expect_equal(nrow(byTarget), 2L)
    byMech2 <- .dtiParseDrugTargets(drg, "CHEMBL25", "mechanism")
    expect_equal(nrow(byMech2), 1L)
    expect_identical(byMech2$approved_symbol, "PTGS1; PTGS2")
    expect_identical(.dtiParseDrugTargets(NULL, "x", "target"), .dtiEmptyTargets())
})

## ------------------------------------------------------------------
## fields = "core"/"all"/<character vector> shared infrastructure
## ------------------------------------------------------------------

test_that(".dtiSelectFields narrows to core/all/custom columns, no network", {
    wide <- data.frame(QueryIDs = "Q1", chembl_id = "C1", Drug_Name = "ASPIRIN",
                       extra.foo = "x", extra.bar = "y", stringsAsFactors = FALSE)
    core <- c("QueryIDs", "chembl_id", "Drug_Name")

    expect_identical(.dtiSelectFields(wide, "core", core), wide[, core])
    expect_identical(.dtiSelectFields(wide, "all", core), wide)

    sub <- .dtiSelectFields(wide, c("Drug_Name", "extra.foo"), core)
    expect_identical(names(sub), c("QueryIDs", "Drug_Name", "extra.foo"))

    expect_error(.dtiSelectFields(wide, "NOT_A_FIELD", core),
                "Unknown field")
})

test_that(".dtiFlattenRecord flattens scalars, nested objects and arrays, no network", {
    rec <- list(
        pref_name = "ASPIRIN",
        max_phase = 4,
        molecule_properties = list(alogp = 1.2, hba = 3),
        target_components = list(
            list(accession = "P00519", component_type = "PROTEIN"),
            list(accession = "P00520", component_type = "PROTEIN"))
    )
    flat <- .dtiFlattenRecord(rec, "molecule")
    expect_identical(flat[["molecule.pref_name"]], "ASPIRIN")
    expect_identical(flat[["molecule.max_phase"]], "4")
    expect_identical(flat[["molecule.molecule_properties.alogp"]], "1.2")
    expect_identical(flat[["molecule.target_components.accession"]],
                     "P00519; P00520")
    expect_identical(.dtiFlattenRecord(NULL, "x"), list())
    expect_identical(.dtiFlattenRecord(list(), "x"), list())
})

test_that(".dtiFlattenGrouped collapses scalar fields across several records, no network", {
    recs <- list(list(mesh_id = "D1", mesh_heading = "Pain"),
                list(mesh_id = "D2", mesh_heading = "Fever"))
    flat <- .dtiFlattenGrouped(recs, "indication")
    expect_identical(flat[["indication.mesh_id"]], "D1; D2")
    expect_identical(flat[["indication.mesh_heading"]], "Pain; Fever")
    expect_identical(.dtiFlattenGrouped(list(), "x"), list())
})

test_that("listDrugTargetFields returns the documented static column list, no network", {
    markers <- list(chembl = "QueryIDs", pubchem = "gene_symbol",
                    dgidb = "gene_name", opentargets = "ensembl_id")
    for (src in names(markers)) {
        f <- listDrugTargetFields(src)
        expect_type(f, "character")
        expect_true(length(f) > 10L)
        expect_true(markers[[src]] %in% f)
    }
    expect_error(listDrugTargetFields("not_a_source"))
})

## ------------------------------------------------------------------
## fields = "all" on the live get*DrugTarget() source functions
## ------------------------------------------------------------------

## The 4 tests below each make exactly one live call (fields = "all")
## rather than a separate live "core" baseline + a separate live custom-
## vector call: fields = "core" is by construction a fixed subset of
## fields = "all" (.dtiSelectFields()), already documented with zero
## network calls via listDrugTargetFields()/listBioassayFields(), and
## already exercised live by each function's own dedicated core test
## elsewhere in this file - re-fetching it here would just be a second
## live round trip to the same endpoint. Custom-vector selection is
## client-side logic on the already-fetched wide data, so it's tested
## directly via .dtiSelectFields() on the "all" result already in hand.

test_that("getChemblDrugTarget(fields = 'all') adds source-prefixed columns", {
    skip_if_offline_dti()
    all <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL1421"), fields = "all")
    coreCols <- c("QueryIDs", "chembl_id", "Drug_Name", "PubChem_CID", "MOA",
                 "Action_Type", "Max_Phase", "First_Approval", "ChEMBL_TID",
                 "UniProt_ID", "Desc", "Organism", "Mesh_Indication")
    expect_identical(coreCols, listDrugTargetFields("chembl")[seq_along(coreCols)])
    expect_true(all(coreCols %in% names(all)))
    expect_true(ncol(all) > length(coreCols))
    expect_true("molecule.molecule_type" %in% names(all))

    sub <- .dtiSelectFields(all, c("Drug_Name", "mechanism.mechanism_comment"), coreCols)
    expect_identical(names(sub), c("QueryIDs", "Drug_Name",
                                   "mechanism.mechanism_comment"))
})

test_that("getPubchemDrugTarget(fields = 'all') adds activity.*-prefixed columns", {
    skip_if_offline_dti()
    all <- getPubchemDrugs("FGFR1", maxCids = 5L, fields = "all")
    coreCols <- c("gene_symbol", "geneid", "target_accession", "cid", "drug_name",
                 "canonical_smiles", "activity_name", "activity_value_uM",
                 "assay_name", "db")
    expect_identical(coreCols, listDrugTargetFields("pubchem")[seq_along(coreCols)])
    expect_true(all(coreCols %in% names(all)))
    expect_true(ncol(all) > length(coreCols))
    expect_true(any(grepl("^activity\\.", names(all))))
})

test_that("getDgidbDrugTarget(fields = 'all') adds drug./gene./source.-prefixed columns", {
    skip_if_offline_dti()
    all <- getDgidbTargets("imatinib", fields = "all")
    coreCols <- c("gene_name", "drug_name", "drug_concept_id", "drug_aliases",
                 "drug_approved", "interaction_types", "directionality",
                 "interaction_score", "evidence_score", "sources", "db")
    expect_identical(coreCols, listDrugTargetFields("dgidb")[seq_along(coreCols)])
    expect_true(all(coreCols %in% names(all)))
    expect_true(ncol(all) > length(coreCols))
    expect_true("drug.id" %in% names(all))
    expect_true("source.citation" %in% names(all))
})

test_that("getOpenTargetsDrugTarget(fields = 'all') adds target./drug.-prefixed columns", {
    skip_if_offline_dti()
    all <- getOpenTargetsTargets("aspirin", expand = "target", fields = "all")
    ## getOpenTargetsTargets() is the drug -> target direction, whose core
    ## columns are .dtiTargetCols - a different slice of
    ## listDrugTargetFields("opentargets") than the target -> drug
    ## direction's .dtiDrugCols, so hardcoded here rather than sliced by
    ## position from the combined list.
    coreCols <- c("chembl_id", "drug_name", "drug_type", "max_clinical_stage",
                 "mechanism_of_action", "action_type", "moa_target_name",
                 "target_id", "approved_symbol")
    expect_true(all(coreCols %in% listDrugTargetFields("opentargets")))
    expect_true(all(coreCols %in% names(all)))
    expect_true(ncol(all) > length(coreCols))
    expect_true("target.biotype" %in% names(all))
    expect_true("drug.description" %in% names(all))
})
