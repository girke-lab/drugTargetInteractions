## Network-free tests for the shared-join-key layer. Every test passes its
## own synthetic HGNC table, so nothing here downloads or queries anything.

.commonIdsHgnc <- function() {
    df <- data.frame(
        hgnc_id = c("HGNC:3688", "HGNC:76", "HGNC:6121", "HGNC:9999"),
        symbol  = c("FGFR1", "ABL1", "KLB", "TWOACC"),
        ensembl_gene_id = c("ENSG00000077782", "ENSG00000097007",
                            "ENSG00000134962", NA_character_),
        stringsAsFactors = FALSE)
    df$prev_symbol  <- list("FLT2", "ABL", character(0), character(0))
    df$alias_symbol <- list(character(0), "JTK7", "KLOTHOB", character(0))
    ## P11362 is FGFR1's; SHAREDACC deliberately belongs to two genes, and
    ## TWOACC deliberately carries two accessions of its own.
    df$uniprot_ids  <- list("P11362", c("P00519", "SHAREDACC"),
                            "SHAREDACC", c("Q11111", "Q22222"))
    df
}

test_that("addCommonIds fills hgnc_id from an accession, an Ensembl id or a symbol", {
    hgnc <- .commonIdsHgnc()
    res <- list(
        chembl      = data.frame(UniProt_ID = c("P11362", "P00519"),
                                 chembl_id = c("CHEMBL941", "CHEMBL1421"),
                                 stringsAsFactors = FALSE),
        opentargets = data.frame(ensembl_id = "ENSG00000077782",
                                 approved_symbol = "FGFR1",
                                 drug_id = "CHEMBL403989", stringsAsFactors = FALSE),
        broad       = data.frame(target_gene = "FLT2", stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))

    expect_identical(out$chembl$hgnc_id, c("HGNC:3688", "HGNC:76"))
    expect_identical(out$chembl$target_uniprot, c("P11362", "P00519"))
    expect_identical(out$chembl$compound_chembl_id, c("CHEMBL941", "CHEMBL1421"))
    expect_identical(out$opentargets$hgnc_id, "HGNC:3688")
    expect_identical(out$opentargets$compound_chembl_id, "CHEMBL403989")
    ## "FLT2" is FGFR1's retired symbol - normalized before lookup.
    expect_identical(out$broad$hgnc_id, "HGNC:3688")
    expect_identical(out$broad$gene_symbol, "FGFR1")
    expect_true(is.na(out$broad$compound_chembl_id))
})

test_that("addCommonIds leaves ambiguous identifiers NA instead of guessing", {
    hgnc <- .commonIdsHgnc()
    res <- list(
        ## SHAREDACC belongs to both ABL1 and KLB - unresolvable.
        chembl = data.frame(UniProt_ID = "SHAREDACC", chembl_id = "CHEMBL1",
                            stringsAsFactors = FALSE),
        ## TWOACC has two accessions, so no single one can be reported.
        broad  = data.frame(target_gene = "TWOACC", stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))

    expect_true(is.na(out$chembl$hgnc_id))
    expect_true(is.na(out$chembl$gene_symbol))
    expect_identical(out$broad$hgnc_id, "HGNC:9999")
    expect_true(is.na(out$broad$target_uniprot))
})

test_that("addCommonIds recognises a ChEMBL id only where the source really has one", {
    hgnc <- .commonIdsHgnc()
    res <- list(dgidb = data.frame(
        gene_name = rep("FGFR1", 3),
        drug_concept_id = c("chembl:CHEMBL52885", "ncit:C104267", "rxcui:357977"),
        stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))
    expect_identical(out$dgidb$compound_chembl_id,
                     c("CHEMBL52885", NA_character_, NA_character_))
})

test_that("addCommonIds never changes the row count or row order", {
    hgnc <- .commonIdsHgnc()
    ## TWOACC's two accessions are exactly the case a merge() would fan out.
    res <- list(broad = data.frame(target_gene = c("TWOACC", "FGFR1", "TWOACC"),
                                   pert_iname = c("a", "b", "c"),
                                   stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))
    expect_identical(nrow(out$broad), 3L)
    expect_identical(out$broad$pert_iname, c("a", "b", "c"))
})

test_that("addCommonIds validates its input and tolerates empty tables", {
    hgnc <- .commonIdsHgnc()
    expect_error(addCommonIds(data.frame(x = 1)), "named list")
    expect_error(addCommonIds(list(notasource = data.frame(x = 1))),
                 "does not recognise source")
    expect_identical(addCommonIds(list()), list())
    empty <- list(chembl = data.frame(UniProt_ID = character(0),
                                      chembl_id = character(0)))
    expect_identical(suppressWarnings(addCommonIds(empty, hgncTable = hgnc)), empty)
})

test_that("combineDrugTargets can carry the shared keys into the appended table", {
    hgnc <- .commonIdsHgnc()
    res <- list(
        chembl = data.frame(QueryIDs = "P11362", UniProt_ID = "P11362",
                            chembl_id = "CHEMBL941", Drug_Name = "IMATINIB",
                            Action_Type = "INHIBITOR", stringsAsFactors = FALSE),
        broad  = data.frame(QueryIDs = "FGFR1", target_gene = "FGFR1",
                            pert_iname = "azd4547", moa = "FGFR inhibitor",
                            stringsAsFactors = FALSE))
    ## Keys pre-computed, so this stays network-free.
    res <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))
    out <- combineDrugTargets(res, columns = c("hgnc_id", "compound_chembl_id",
                                               "drug_name", "source"))
    expect_identical(names(out), c("hgnc_id", "compound_chembl_id",
                                   "drug_name", "source"))
    expect_identical(out$hgnc_id, c("HGNC:3688", "HGNC:3688"))
    expect_identical(out$compound_chembl_id, c("CHEMBL941", NA_character_))
    expect_identical(out$source, c("ChEMBL", "Broad Repurposing Hub"))

    ## ChEMBL has no gene-symbol column of its own, so gene_symbol is NA for
    ## its rows by default - but once the shared keys are present it comes
    ## from there, with no network call.
    withKeys <- combineDrugTargets(res, columns = c("hgnc_id", "gene_symbol", "source"))
    expect_identical(withKeys$gene_symbol, c("FGFR1", "FGFR1"))
    plain <- combineDrugTargets(res["broad"], columns = c("gene_symbol", "source"))
    expect_identical(plain$gene_symbol, "FGFR1")
})
