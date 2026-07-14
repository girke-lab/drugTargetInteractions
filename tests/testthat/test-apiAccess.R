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
