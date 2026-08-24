## Network-guarded tests for the HGNC-anchored genome-wide master table
## builder. Mirrors the no-mocking, skip-when-offline convention used
## throughout this package's test suite.

## skip_on_bioc() is a deliberate, temporary over-correction for this
## major-upgrade release - see test-apiAccess.R's copy of this helper
## for the full rationale. Remove selectively later, not all at once.
skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## --- getHgncGeneTable() / downloadHgncTable() -------------------------

test_that("getHgncGeneTable returns a clean, gene-centric protein-coding table", {
    skip_if_offline_dti()
    hgncTable <- getHgncGeneTable()
    expect_s3_class(hgncTable, "data.frame")
    expect_true(nrow(hgncTable) > 15000L)  # ~19,200 as of 2026-07
    expect_identical(names(hgncTable),
                     c("hgnc_id", "symbol", "prev_symbol", "alias_symbol",
                       "entrez_id", "ensembl_gene_id", "uniprot_ids", "locus_group"))
    expect_true(all(hgncTable$locus_group == "protein-coding gene"))
    expect_true(is.list(hgncTable$uniprot_ids))
    expect_true(is.list(hgncTable$prev_symbol))
    expect_true(is.list(hgncTable$alias_symbol))

    fgfr1 <- hgncTable[hgncTable$symbol == "FGFR1", ]
    expect_equal(nrow(fgfr1), 1L)
    expect_identical(fgfr1$hgnc_id, "HGNC:3688")
    expect_identical(fgfr1$ensembl_gene_id, "ENSG00000077782")
    expect_identical(fgfr1$uniprot_ids[[1]], "P11362")

    ## No literal quote characters should leak into any multi-value
    ## column (regression check for the CSV-style-quoting bug fixed this
    ## session: read.delim(quote = "") left wrapping quotes in the data
    ## instead of stripping them).
    expect_false(any(grepl("\"", unlist(hgncTable$alias_symbol), fixed = TRUE)))
    expect_false(any(grepl("\"", unlist(hgncTable$prev_symbol), fixed = TRUE)))
    expect_false(any(grepl("\"", unlist(hgncTable$uniprot_ids), fixed = TRUE)))
})

test_that("getHgncGeneTable(proteinCodingOnly = FALSE) keeps non-coding loci too", {
    skip_if_offline_dti()
    full <- getHgncGeneTable(proteinCodingOnly = FALSE)
    coding <- getHgncGeneTable(proteinCodingOnly = TRUE)
    expect_true(nrow(full) > nrow(coding))
})

## --- buildHgncSymbolMap() / normalizeGeneSymbols() (no network) -------

.syntheticHgncTable <- function() {
    data.frame(
        symbol = c("FGFR1", "ABL1", "KLB"),
        stringsAsFactors = FALSE
    ) -> df
    df$prev_symbol <- list("FLT2", "ABL", character(0))
    df$alias_symbol <- list(character(0), "JTK7", c("SHARED", "KLOTHOB"))
    ## Deliberately ambiguous: "SHARED" is also an alias of ABL1, so it
    ## maps to two different current symbols.
    df$alias_symbol[[2]] <- c("JTK7", "SHARED")
    df
}

test_that("buildHgncSymbolMap resolves old symbols and flags ambiguity, no network", {
    hgncTable <- .syntheticHgncTable()
    expect_warning(map <- buildHgncSymbolMap(hgncTable), "ambiguous")
    expect_identical(unname(map["FLT2"]), "FGFR1")
    expect_identical(unname(map["ABL"]), "ABL1")
    expect_identical(unname(map["JTK7"]), "ABL1")
    expect_identical(unname(map["KLOTHOB"]), "KLB")
    ## "SHARED" is ambiguous (ABL1 and KLB both claim it) - resolved
    ## deterministically to the alphabetically-first current symbol.
    expect_identical(unname(map["SHARED"]), "ABL1")
    expect_true("SHARED" %in% names(attr(map, "ambiguous")))
    expect_setequal(attr(map, "ambiguous")[["SHARED"]], c("ABL1", "KLB"))
})

test_that("buildHgncSymbolMap handles a table with no prev/alias symbols at all", {
    hgncTable <- data.frame(symbol = "FGFR1", stringsAsFactors = FALSE)
    hgncTable$prev_symbol <- list(character(0))
    hgncTable$alias_symbol <- list(character(0))
    map <- buildHgncSymbolMap(hgncTable)
    expect_length(map, 0L)
    expect_length(attr(map, "ambiguous"), 0L)
})

test_that("buildHgncSymbolMap can be built without the map-wide warning", {
    hgncTable <- .syntheticHgncTable()
    ## The count describes the HGNC snapshot, so callers translating a
    ## handful of symbols suppress it and report what those symbols hit.
    expect_silent(map <- buildHgncSymbolMap(hgncTable, warn = FALSE))
    ## Suppressing the warning must not change the map or its attribute.
    expect_identical(map, suppressWarnings(buildHgncSymbolMap(hgncTable)))
    expect_true("SHARED" %in% names(attr(map, "ambiguous")))
})

test_that("normalizeGeneSymbols passes through current symbols and resolves old ones", {
    hgncTable <- .syntheticHgncTable()
    ## None of these is ambiguous, so translating them says nothing - the
    ## map holding an ambiguous entry elsewhere is not this caller's problem.
    expect_silent(out <- normalizeGeneSymbols(
        c("FGFR1", "FLT2", "ABL", "NOTAREALSYMBOL"), hgncTable = hgncTable))
    expect_identical(as.character(out), c("FGFR1", "FGFR1", "ABL1", NA_character_))
    expect_identical(attr(out, "unmapped"), "NOTAREALSYMBOL")
    expect_length(attr(out, "ambiguous"), 0L)
})

test_that("normalizeGeneSymbols warns only about the symbols it was given", {
    hgncTable <- .syntheticHgncTable()
    ## "SHARED" is claimed by both ABL1 and KLB.
    expect_warning(out <- normalizeGeneSymbols(c("FGFR1", "SHARED"),
                                               hgncTable = hgncTable),
                   "1 of the symbol\\(s\\) given")
    expect_identical(as.character(out), c("FGFR1", "ABL1"))
    ## The reported ambiguity is scoped to what was asked about, and names it.
    expect_identical(names(attr(out, "ambiguous")), "SHARED")
    expect_setequal(attr(out, "ambiguous")[["SHARED"]], c("ABL1", "KLB"))
})

test_that("a symbol that is already current is never reported as ambiguous", {
    hgncTable <- .syntheticHgncTable()
    ## Make "SHARED" both an approved symbol and an ambiguous alias: it is
    ## passed through untouched, so the alias collision cannot apply to it.
    hgncTable[4, "symbol"] <- "SHARED"
    hgncTable$prev_symbol[[4]] <- character(0)
    hgncTable$alias_symbol[[4]] <- character(0)
    expect_silent(out <- normalizeGeneSymbols("SHARED", hgncTable = hgncTable))
    expect_identical(as.character(out), "SHARED")
})

## --- buildGenomeWideDrugTargetTable() checkpointing/resume ------------

test_that("buildGenomeWideDrugTargetTable rejects a missing outDir without any network calls", {
    hgncTable <- .syntheticHgncTable()
    hgncTable$ensembl_gene_id <- c("ENSG1", "ENSG2", "ENSG3")
    hgncTable$hgnc_id <- c("HGNC:1", "HGNC:2", "HGNC:3")
    hgncTable$uniprot_ids <- list("P11362", "P00519", "Q86Z14")
    expect_error(
        buildGenomeWideDrugTargetTable(hgncTable = hgncTable, sources = "dgidb"),
        "outDir")
    expect_error(
        buildGenomeWideDrugTargetTable(hgncTable = hgncTable, sources = "ttd",
                                       outDir = tempfile()),
        "ttdDbPath")
})

test_that("buildGenomeWideDrugTargetTable checkpoints per chunk and resumes correctly", {
    skip_if_offline_dti()
    hgncTable <- getHgncGeneTable()
    hgncSmall <- hgncTable[hgncTable$symbol %in% c("FGFR1", "KLB", "TFEB"), ]
    outDir <- tempfile("dti_build_test_")
    on.exit({
        bfc <- .getCache()
        r <- bfcquery(bfc, "genomewide-drugtargets-", field = "rname")
        if (nrow(r) > 0) bfcremove(bfc, r$rid)
        unlink(outDir, recursive = TRUE)
    }, add = TRUE)

    res1 <- buildGenomeWideDrugTargetTable(
        hgncTable = hgncSmall, sources = "dgidb",
        outDir = outDir, chunkGenes = 1L, verbose = FALSE)
    expect_true(nrow(res1$dgidb) > 0L)
    expect_true(all(c("hgnc_id", "symbol", "ensembl_gene_id", "QueryIDs") %in%
                    names(res1$dgidb)))
    manifest <- readRDS(file.path(outDir, "manifest.rds"))
    expect_length(manifest, 3L)  # one chunk per gene, chunkGenes = 1

    ## Simulate an interrupted chunk, then confirm resuming only re-runs
    ## the missing one and reproduces the identical final result.
    file.remove(file.path(outDir, "dgidb_chunk2.rds"))
    saveRDS(setdiff(manifest, "dgidb_chunk2"), file.path(outDir, "manifest.rds"))
    res2 <- buildGenomeWideDrugTargetTable(
        hgncTable = hgncSmall, sources = "dgidb",
        outDir = outDir, chunkGenes = 1L, verbose = FALSE, rerun = FALSE)
    expect_identical(res1$dgidb, res2$dgidb)

    ## The cache name must not collide with (be a substring match for)
    ## an unrelated rname - see the cacheName construction in
    ## buildGenomeWideDrugTargetTable() for why this matters (bfcquery()'s
    ## default rname matching is substring-based, not exact).
    expect_false(is.null(attr(res2, "cachePath")))
})
