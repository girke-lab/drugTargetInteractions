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

## --- mergeDrugTargets() (horizontal / column append) -------------------

.mergeFixture <- function() {
    hgnc <- .commonIdsHgnc()
    res <- list(
        ## two mechanisms for one drug, plus a second drug - the many-per-key
        ## case the collapse exists for
        chembl = data.frame(
            QueryIDs = "P11362", UniProt_ID = "P11362",
            chembl_id = c("CHEMBL941", "CHEMBL941", "CHEMBL1421"),
            Drug_Name = c("IMATINIB", "IMATINIB", "DASATINIB"),
            Action_Type = c("INHIBITOR", "BLOCKER", "INHIBITOR"),
            stringsAsFactors = FALSE),
        opentargets = data.frame(
            QueryIDs = "FGFR1", ensembl_id = "ENSG00000077782",
            approved_symbol = "FGFR1", drug_id = "CHEMBL941",
            drug_name = "IMATINIB", action_type = "INHIBITOR",
            stringsAsFactors = FALSE),
        ## no compound id of its own - contributes only under by = "hgnc_id"
        broad = data.frame(QueryIDs = "FGFR1", target_gene = "FGFR1",
                           pert_iname = c("azd4547", "brivanib"),
                           stringsAsFactors = FALSE))
    suppressWarnings(addCommonIds(res, hgncTable = hgnc))
}

test_that("mergeDrugTargets puts sources side by side, one row per gene-drug pair", {
    out <- mergeDrugTargets(.mergeFixture(),
                            columns = c("Drug_Name", "Action_Type", "drug_name"))
    expect_s4_class(out, "DataFrame")
    ## CHEMBL941 and CHEMBL1421, both for FGFR1 - Broad has no compound id
    expect_identical(nrow(out), 2L)
    expect_identical(out$compound_chembl_id, c("CHEMBL1421", "CHEMBL941"))
    imatinib <- which(out$compound_chembl_id == "CHEMBL941")
    expect_identical(out$n_sources[imatinib], 2L)
    expect_setequal(out$sources[[imatinib]], c("chembl", "opentargets"))
    ## the two ChEMBL mechanisms collapsed into one cell, drug name deduplicated
    expect_setequal(out$chembl_Action_Type[[imatinib]], c("INHIBITOR", "BLOCKER"))
    expect_identical(out$chembl_Drug_Name[[imatinib]], "IMATINIB")
    ## a source absent from a key gets an empty cell, not a missing column
    dasatinib <- which(out$compound_chembl_id == "CHEMBL1421")
    expect_identical(out$opentargets_drug_name[[dasatinib]], character(0))
})

test_that("mergeDrugTargets keyed on the gene keeps sources that have no compound id", {
    out <- mergeDrugTargets(.mergeFixture(), by = "hgnc_id",
                            columns = c("Drug_Name", "pert_iname"))
    expect_identical(nrow(out), 1L)
    expect_identical(out$n_sources, 3L)
    expect_setequal(out$chembl_Drug_Name[[1]], c("IMATINIB", "DASATINIB"))
    expect_setequal(out$broad_pert_iname[[1]], c("azd4547", "brivanib"))
})

test_that("mergeDrugTargets reports the rows it drops for want of a key", {
    expect_message(
        mergeDrugTargets(.mergeFixture(), columns = "Drug_Name", verbose = TRUE),
        "broad - dropped 2/2 row\\(s\\) with no hgnc_id/compound_chembl_id")
})

test_that("mergeDrugTargets collapse = 'string' pastes each cell", {
    out <- mergeDrugTargets(.mergeFixture(), by = "hgnc_id",
                            columns = "Action_Type", collapse = "string",
                            sep = " | ")
    expect_type(out$chembl_Action_Type, "character")
    expect_identical(out$chembl_Action_Type, "INHIBITOR | BLOCKER")
})

test_that("mergeDrugTargets validates input and survives having nothing to join", {
    expect_error(mergeDrugTargets(data.frame(x = 1)), "named list")
    expect_error(mergeDrugTargets(list(notasource = data.frame(x = 1))),
                 "does not recognise source")
    ## every row lacking a key leaves an empty table with the key columns
    onlyBroad <- .mergeFixture()["broad"]
    out <- mergeDrugTargets(onlyBroad)
    expect_identical(nrow(out), 0L)
    expect_identical(names(out), c("hgnc_id", "compound_chembl_id"))
})

test_that("addCommonIds keeps an hgnc_id it was given and only fills the gaps", {
    hgnc <- .commonIdsHgnc()
    ## Row 1 arrives already keyed, row 2 does not. Re-deriving row 1 from
    ## its symbol would move it to FGFR1; a genome-wide build records the
    ## gene each row was queried from, so what it supplies must stand.
    res <- list(ttd = data.frame(
        hgnc_id = c("HGNC:6121", NA_character_),
        GeneName = c("FGFR1", "ABL1"), Uniprot_acc = NA_character_,
        stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))

    expect_identical(out$ttd$hgnc_id, c("HGNC:6121", "HGNC:76"))
    ## The dependent columns follow the hgnc_id that was kept.
    expect_identical(out$ttd$gene_symbol, c("KLB", "ABL1"))
    expect_identical(nrow(out$ttd), 2L)
})

test_that("addCommonIds is unchanged for tables that carry no hgnc_id", {
    hgnc <- .commonIdsHgnc()
    res <- list(chembl = data.frame(UniProt_ID = c("P11362", "P00519"),
                                    chembl_id = c("CHEMBL941", "CHEMBL1421"),
                                    stringsAsFactors = FALSE))
    out <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))
    expect_identical(out$chembl$hgnc_id, c("HGNC:3688", "HGNC:76"))
})

test_that("addCommonIds is idempotent", {
    hgnc <- .commonIdsHgnc()
    res <- list(chembl = data.frame(UniProt_ID = "P11362", chembl_id = "CHEMBL941",
                                    stringsAsFactors = FALSE))
    once  <- suppressWarnings(addCommonIds(res, hgncTable = hgnc))
    twice <- suppressWarnings(addCommonIds(once, hgncTable = hgnc))
    expect_identical(once, twice)
})

## mergeDrugTargets() derives keys with addCommonIds() only when a table
## lacks them, and addCommonIds() takes its HGNC table as an argument
## while mergeDrugTargets() does not. The tests below therefore hand it
## tables that already carry all four key columns, which keeps them
## network-free - and is what a genome-wide build supplies anyway.

test_that("mergeDrugTargets carries gene identity once, not once per source", {
    ## Both sources tag their rows with the same gene identity, the way a
    ## genome-wide build does.
    ident <- function(...) data.frame(
        hgnc_id = "HGNC:3688", symbol = "FGFR1",
        ensembl_gene_id = "ENSG00000077782", gene_symbol = "FGFR1",
        target_uniprot = "P11362", ..., stringsAsFactors = FALSE)
    res <- list(
        chembl = ident(compound_chembl_id = "CHEMBL941", Drug_Name = "IMATINIB"),
        ttd    = ident(compound_chembl_id = NA_character_, DrugName = "Debio 1347"))
    out <- mergeDrugTargets(res, by = "hgnc_id")

    expect_identical(names(out)[1:5],
                     c("hgnc_id", "symbol", "ensembl_gene_id", "n_sources", "sources"))
    expect_identical(out$symbol[[1]], "FGFR1")
    expect_identical(out$ensembl_gene_id[[1]], "ENSG00000077782")
    ## No per-source copies of what the key already determines.
    expect_length(grep("_(symbol|ensembl_gene_id)$", names(out)), 0L)
    ## The values themselves still arrive, per source.
    expect_identical(out$chembl_Drug_Name[[1]], "IMATINIB")
    expect_identical(out$ttd_DrugName[[1]], "Debio 1347")
})

test_that("mergeDrugTargets keeps every gene when the key is a compound", {
    ## One drug against two genes: keyed on the compound, `symbol` is
    ## legitimately both of them rather than one arbitrary winner.
    res <- list(chembl = data.frame(
        hgnc_id = c("HGNC:3688", "HGNC:76"), symbol = c("FGFR1", "ABL1"),
        ensembl_gene_id = c("ENSG00000077782", "ENSG00000097007"),
        gene_symbol = c("FGFR1", "ABL1"), target_uniprot = c("P11362", "P00519"),
        compound_chembl_id = "CHEMBL941", Drug_Name = "IMATINIB",
        stringsAsFactors = FALSE))
    out <- mergeDrugTargets(res, by = "compound_chembl_id")

    expect_identical(nrow(out), 1L)
    expect_identical(sort(out$symbol[[1]]), c("ABL1", "FGFR1"))
    expect_identical(sort(out$ensembl_gene_id[[1]]),
                     c("ENSG00000077782", "ENSG00000097007"))
})

test_that("mergeDrugTargets resolves a shared column name to each source's own", {
    keyed <- function(...) data.frame(
        hgnc_id = "HGNC:3688", gene_symbol = "FGFR1", target_uniprot = "P11362",
        ..., stringsAsFactors = FALSE)
    res <- list(
        chembl = keyed(compound_chembl_id = "CHEMBL941", Drug_Name = "IMATINIB",
                       Action_Type = "INHIBITOR", Max_Phase = 4L),
        ttd    = keyed(compound_chembl_id = NA_character_, DrugName = "Debio 1347",
                       MOA = "Inhibitor", Highest_status = "Phase 2"))

    ## "drug_name" picks up Drug_Name and DrugName without naming either.
    byShared <- mergeDrugTargets(res, by = "hgnc_id", columns = "drug_name")
    expect_true(all(c("chembl_Drug_Name", "ttd_DrugName") %in% names(byShared)))
    expect_false(any(grepl("Action_Type|Max_Phase", names(byShared))))

    ## Naming each source's own column still works exactly as before.
    byNative <- mergeDrugTargets(res, by = "hgnc_id",
                                 columns = c("Drug_Name", "DrugName"))
    expect_identical(names(byShared), names(byNative))

    ## The two forms mix freely, and a shared name spanning two concepts
    ## brings each source's variant of both.
    both <- mergeDrugTargets(res, by = "hgnc_id",
                             columns = c("drug_name", "action", "Max_Phase"))
    expect_true(all(c("chembl_Drug_Name", "ttd_DrugName", "chembl_Action_Type",
                      "ttd_MOA", "chembl_Max_Phase") %in% names(both)))
})
