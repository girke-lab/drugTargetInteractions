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
