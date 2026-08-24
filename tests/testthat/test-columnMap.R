## Network-free tests for the column mapping table and the row-append
## assembler that consumes it. Nothing here downloads or queries anything:
## the tests that need shared keys pass their own synthetic HGNC table,
## and the rest run with keys = FALSE.

.colMapHgnc <- function() {
    df <- data.frame(
        hgnc_id = c("HGNC:3688", "HGNC:76", "HGNC:6121"),
        symbol  = c("FGFR1", "ABL1", "KLB"),
        ensembl_gene_id = c("ENSG00000077782", "ENSG00000097007",
                            "ENSG00000134962"),
        stringsAsFactors = FALSE)
    df$prev_symbol  <- list("FLT2", "ABL", character(0))
    df$alias_symbol <- list(character(0), "JTK7", "KLOTHOB")
    df$uniprot_ids  <- list("P11362", "P00519", "Q86Z14")
    df
}

## A two-source genome-wide-shaped build: the identity block a real build
## attaches, plus each source's own columns.
.colMapBuild <- function() {
    list(
        chembl = data.frame(
            QueryIDs = c("P11362", "P11362"), hgnc_id = "HGNC:3688",
            symbol = "FGFR1", ensembl_gene_id = "ENSG00000077782",
            chembl_id = c("CHEMBL941", "CHEMBL1421"),
            Drug_Name = c("IMATINIB", "ERDAFITINIB"),
            Action_Type = c("INHIBITOR", "INHIBITOR"),
            MOA = c("Abl inhibitor", "FGFR inhibitor"),
            Max_Phase = c(4L, 4L), UniProt_ID = "P11362",
            ChEMBL_TID = c("CHEMBL1862", "CHEMBL3650"),
            stringsAsFactors = FALSE),
        ttd = data.frame(
            QueryIDs = "ABL1", hgnc_id = "HGNC:76", symbol = "ABL1",
            ensembl_gene_id = "ENSG00000097007",
            GeneName = "ABL1", DrugName = "Debio 1347", MOA = "Inhibitor",
            Highest_status = "Phase 2", Uniprot_acc = "P00519",
            Smiles = "CC(C)c1ccc", stringsAsFactors = FALSE))
}


## ---------------------------------------------------------------------
## The map itself
## ---------------------------------------------------------------------

test_that("drugTargetColumnMap returns a well-formed, all-active curated map", {
    m <- drugTargetColumnMap()
    expect_s3_class(m, "data.frame")
    expect_identical(names(m), c("canonical", "source", "column", "status", "active"))
    expect_true(all(m$status == "curated"))
    expect_true(all(m$active))
    expect_true(all(m$source %in% c("chembl", "dgidb", "opentargets", "ttd",
                                    "broad", "gtopdb")))
    ## One row per (canonical, source) - anything else is ambiguous.
    expect_false(any(duplicated(paste(m$canonical, m$source))))
    ## The keys addCommonIds() owns are deliberately not mapped.
    expect_false(any(c("hgnc_id", "compound_chembl_id") %in% m$canonical))
})

test_that("the map the user sees is the one combineDrugTargets applies", {
    ## .dtiCombineColMap is derived from the curated map, so the mapping
    ## shown by drugTargetColumnMap() cannot drift from the one actually
    ## used by the older row-append function.
    m <- drugTargetColumnMap()
    legacy <- drugTargetInteractions:::.dtiCombineColMap
    for (canon in names(legacy)) {
        rows <- m[m$canonical == canon, ]
        expect_identical(stats::setNames(rows$column, rows$source), legacy[[canon]])
    }
})

test_that("a native column may feed more than one canonical column", {
    m <- drugTargetColumnMap()
    ## The Broad Hub's free-text `moa` is both the coarse action label and
    ## the finer mechanism descriptor.
    expect_identical(sort(m$canonical[m$source == "broad" & m$column == "moa"]),
                     c("action", "mechanism"))
})


## ---------------------------------------------------------------------
## Proposals
## ---------------------------------------------------------------------

test_that("name-matched groups are proposed but never active", {
    res <- list(
        opentargets = data.frame(approved_symbol = "FGFR1", drug_id = "CHEMBL1",
                                 stringsAsFactors = FALSE),
        ttd = data.frame(GeneName = "FGFR1", DrugID = "D03DSN",
                         stringsAsFactors = FALSE))
    m <- drugTargetColumnMap(res)
    prop <- m[m$status == "proposed", ]

    ## drug_id/DrugID normalise to the same name but are different
    ## identifier spaces - the proposal must not take effect on its own.
    expect_identical(sort(prop$column), c("DrugID", "drug_id"))
    expect_false(any(prop$active))
    ## ...and it is genuinely inert until switched on.
    expect_false("drugid" %in% names(combineGenomeWideDrugTargets(res, keys = FALSE)))
    m$active[m$status == "proposed"] <- TRUE
    expect_true("drugid" %in% names(combineGenomeWideDrugTargets(res, colMap = m,
                                                                keys = FALSE)))
})

test_that("a column in no group at all is reported as unmatched", {
    res <- list(ttd = data.frame(GeneName = "FGFR1", Smiles = "CC",
                                 stringsAsFactors = FALSE))
    m <- drugTargetColumnMap(res)
    expect_identical(m$column[m$status == "unmatched"], "Smiles")
    expect_true(is.na(m$canonical[m$status == "unmatched"]))
    ## A name occurring in only one source is not a cross-source group.
    expect_identical(nrow(m[m$status == "proposed", ]), 0L)
})

test_that("proposals skip columns the curated map already names", {
    res <- list(chembl = data.frame(Drug_Name = "IMATINIB", stringsAsFactors = FALSE),
                ttd    = data.frame(DrugName = "IMATINIB", stringsAsFactors = FALSE))
    m <- drugTargetColumnMap(res)
    expect_identical(nrow(m[m$status != "curated", ]), 0L)
})

test_that("the identity block and the derived keys are never proposed", {
    res <- list(chembl = data.frame(hgnc_id = "HGNC:3688", symbol = "FGFR1",
                                    ensembl_gene_id = "ENSG00000077782",
                                    QueryIDs = "P11362",
                                    compound_chembl_id = "CHEMBL941",
                                    stringsAsFactors = FALSE),
                ttd    = data.frame(hgnc_id = "HGNC:76", symbol = "ABL1",
                                    ensembl_gene_id = "ENSG00000097007",
                                    QueryIDs = "ABL1",
                                    compound_chembl_id = NA_character_,
                                    stringsAsFactors = FALSE))
    expect_identical(nrow(drugTargetColumnMap(res)[
        drugTargetColumnMap(res)$status != "curated", ]), 0L)
})


## ---------------------------------------------------------------------
## Validation of an edited map
## ---------------------------------------------------------------------

test_that("an edited map is rejected when it cannot be applied sensibly", {
    res <- .colMapBuild()

    expect_error(combineGenomeWideDrugTargets(res, colMap = "not a map"),
                 "must be a data.frame")
    expect_error(combineGenomeWideDrugTargets(
        res, colMap = data.frame(canonical = "drug_name", source = "ttd")),
        "missing required column")

    ## One source cannot fill one canonical column from two columns.
    m <- drugTargetColumnMap()
    m <- rbind(m, data.frame(canonical = "drug_name", source = "ttd",
                             column = "Smiles", status = "curated", active = TRUE))
    expect_error(combineGenomeWideDrugTargets(res, colMap = m, keys = FALSE),
                 "ambiguous")

    ## The cross-source keys are not the map's to assign.
    m2 <- drugTargetColumnMap()
    m2 <- rbind(m2, data.frame(canonical = "hgnc_id", source = "ttd",
                               column = "GeneName", status = "curated", active = TRUE))
    expect_error(combineGenomeWideDrugTargets(res, colMap = m2, keys = FALSE),
                 "addCommonIds")
})

test_that("deactivating a row drops that column without affecting others", {
    res <- .colMapBuild()
    m <- drugTargetColumnMap()
    m$active[m$canonical == "max_phase"] <- FALSE
    out <- combineGenomeWideDrugTargets(res, colMap = m, keys = FALSE)
    expect_false("max_phase" %in% names(out))
    expect_true(all(c("drug_name", "action") %in% names(out)))
    ## Deactivated, the native columns fall through to the passthrough block.
    expect_true(all(c("chembl.Max_Phase", "ttd.Highest_status") %in% names(out)))
})


## ---------------------------------------------------------------------
## The assembler
## ---------------------------------------------------------------------

test_that("combineGenomeWideDrugTargets appends rows without changing the grain", {
    res <- .colMapBuild()
    out <- combineGenomeWideDrugTargets(res, keys = FALSE)

    expect_identical(nrow(out), sum(vapply(res, nrow, integer(1))))
    expect_identical(as.vector(table(out$source)[c("ChEMBL", "TTD")]), c(2L, 1L))
    ## The identity block comes first and is carried through verbatim.
    expect_identical(names(out)[1:5],
                     c("hgnc_id", "symbol", "ensembl_gene_id", "QueryIDs", "source"))
    expect_identical(out$hgnc_id, c("HGNC:3688", "HGNC:3688", "HGNC:76"))
    expect_identical(out$symbol, c("FGFR1", "FGFR1", "ABL1"))
})

test_that("columns holding the same content under different names are aligned", {
    out <- combineGenomeWideDrugTargets(.colMapBuild(), keys = FALSE)
    expect_identical(out$drug_name, c("IMATINIB", "ERDAFITINIB", "Debio 1347"))
    ## action takes ChEMBL's Action_Type and TTD's MOA...
    expect_identical(out$action, c("INHIBITOR", "INHIBITOR", "Inhibitor"))
    ## ...while mechanism takes ChEMBL's MOA, which TTD has no counterpart for.
    expect_identical(out$mechanism, c("Abl inhibitor", "FGFR inhibitor", NA))
    expect_identical(out$target_uniprot, c("P11362", "P11362", "P00519"))
})

test_that("a shared concept with unshared vocabularies is carried as text", {
    out <- combineGenomeWideDrugTargets(.colMapBuild(), keys = FALSE)
    ## ChEMBL's numeric phase and TTD's phrase land in one column, readable
    ## alongside `source` rather than forced into a common scale.
    expect_type(out$max_phase, "character")
    expect_identical(out$max_phase, c("4", "4", "Phase 2"))
})

test_that("a column belonging to one source is carried through, source-prefixed", {
    out <- combineGenomeWideDrugTargets(.colMapBuild(), keys = FALSE)
    expect_true(all(c("chembl.chembl_id", "chembl.ChEMBL_TID", "ttd.Smiles")
                    %in% names(out)))
    expect_identical(out$ttd.Smiles, c(NA, NA, "CC(C)c1ccc"))
    expect_identical(out$chembl.chembl_id, c("CHEMBL941", "CHEMBL1421", NA))
})

test_that("native = FALSE keeps only the shared columns", {
    out <- combineGenomeWideDrugTargets(.colMapBuild(), native = FALSE, keys = FALSE)
    expect_false(any(grepl("\\.", names(out))))
    expect_true(all(c("drug_name", "action", "max_phase") %in% names(out)))
    expect_identical(nrow(out), 3L)
})

test_that("a source's own column keeps its own type", {
    res <- .colMapBuild()
    res$chembl$First_Approval <- c(2001L, 2019L)
    res$ttd$approved <- TRUE
    out <- combineGenomeWideDrugTargets(res, keys = FALSE)
    ## Each native column comes from exactly one source, so nothing has to
    ## be reconciled and the type survives.
    expect_type(out$chembl.First_Approval, "integer")
    expect_type(out$ttd.approved, "logical")
})

test_that("the gene symbol HGNC gives and the one the source reports are both kept", {
    res <- .colMapBuild()
    res$ttd$GeneName <- "Abl1"          # a source echoing a differently-cased symbol
    out <- combineGenomeWideDrugTargets(res, keys = FALSE)
    expect_identical(out$symbol[3], "ABL1")
    expect_identical(out$gene_symbol[3], "Abl1")
})

test_that("combineGenomeWideDrugTargets rejects a list it cannot interpret", {
    expect_error(combineGenomeWideDrugTargets(data.frame(a = 1)),
                 "named list of per-source")
    expect_error(combineGenomeWideDrugTargets(list(nosuchsource = data.frame(a = 1))),
                 "does not recognise source")
    expect_identical(nrow(combineGenomeWideDrugTargets(list())), 0L)
    expect_identical(nrow(combineGenomeWideDrugTargets(
        list(ttd = data.frame(GeneName = character(0))), keys = FALSE)), 0L)
})

test_that("keys = TRUE derives compound_chembl_id without touching the build's hgnc_id", {
    res <- .colMapBuild()
    ## A deliberately wrong-but-supplied hgnc_id: the assembler must carry
    ## it through rather than silently correcting it from ChEMBL's
    ## accession, so what a build established is what the table reports.
    res$chembl$hgnc_id <- "HGNC:6121"
    out <- suppressWarnings(combineGenomeWideDrugTargets(
        res, hgncTable = .colMapHgnc()))
    expect_identical(out$hgnc_id, c("HGNC:6121", "HGNC:6121", "HGNC:76"))
    expect_identical(out$compound_chembl_id, c("CHEMBL941", "CHEMBL1421", NA))
})
