## Network-guarded tests for the live ChEMBL REST accessors. All tests
## skip cleanly when offline so the Bioconductor build machines never
## fail on a flaky endpoint (BiocCheck penalises unguarded network calls).

skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
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
    expect_identical(names(df), c("QueryIDs", "chembl_id", "Drug_Name", "MOA",
                                  "Action_Type", "Max_Phase", "First_Approval",
                                  "ChEMBL_TID", "UniProt_ID", "Desc",
                                  "Organism", "Mesh_Indication"))
    expect_true(nrow(df) >= 1L)
    expect_true(all(df$QueryIDs == "P11362"))
    expect_true("CHEMBL1201733" %in% df$chembl_id)  # pazopanib
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
                                  "drug_approved", "interaction_types",
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
        drug = list(name = "ASPIRIN", conceptId = "chembl:CHEMBL25", approved = TRUE),
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
    expect_identical(.dtiParseDgidbInteractions(NULL), .dtiEmptyDgidb())
})
