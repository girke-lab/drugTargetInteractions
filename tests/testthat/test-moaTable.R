## Tests for the cross-source MOA assembly layer (assembleMoaTable() and
## its reshapers). assembleMoaTable() is a pure transformer over
## queryDrugTargets()'s output, so the bulk of the coverage here is
## network-free and runs against synthetic per-source frames - that is
## deliberate: it exercises the grain/kind invariants that matter most
## without depending on six live APIs agreeing to be up at once.
##
## The live tests at the end are known-answer checks against real
## upstream content, network-guarded with the shared
## skip_if_offline_dti() pattern - defined locally so this file runs
## standalone (see test-apiAccess.R's copy for the skip_on_bioc()
## rationale).

skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

## --- synthetic per-source fixtures ------------------------------------
## Column names and value shapes mirror what each source's accessor
## really returns (verified live), trimmed to the columns the MOA
## extractor reads.

fxChembl <- function() data.frame(
    QueryIDs = c("CHEMBL941", "CHEMBL941"),
    chembl_id = c("CHEMBL941", "CHEMBL941"),
    Drug_Name = c("IMATINIB", "IMATINIB"),
    MOA = c("Tyrosine-protein kinase ABL inhibitor",
            "Stem cell growth factor receptor inhibitor"),
    Action_Type = c("INHIBITOR", "INHIBITOR"),
    UniProt_ID = c("P00519", "P10721"),
    stringsAsFactors = FALSE)

fxOpenTargets <- function() data.frame(
    QueryIDs = c("imatinib", "imatinib"),
    chembl_id = c("CHEMBL941", "CHEMBL941"),
    drug_name = c("IMATINIB", "IMATINIB"),
    mechanism_of_action = c("Bcr/Abl fusion protein inhibitor",
                            "Bcr/Abl fusion protein inhibitor"),
    action_type = c("INHIBITOR", "INHIBITOR"),
    approved_symbol = c("ABL1", "BCR"),
    stringsAsFactors = FALSE)

## Broad repeats one " | "-joined drug-level string on every target row.
fxBroad <- function() data.frame(
    QueryIDs = c("imatinib", "imatinib", "imatinib"),
    pert_iname = rep("imatinib", 3),
    InChIKey = rep("KTUFNOKKBVMGRW-UHFFFAOYSA-N", 3),
    moa = rep("Bcr-Abl kinase inhibitor | KIT inhibitor", 3),
    target_gene = c("BCR", "ABL1", "CSF1R"),
    stringsAsFactors = FALSE)

fxTtd <- function() data.frame(
    QueryIDs = c("Imatinib", "Imatinib"),
    DrugName = c("Imatinib", "Imatinib"),
    GeneName = c("BCR-ABL1", "MCL1"),
    Uniprot_acc = c("P00519", NA_character_),
    MOA = c("Inhibitor", NA_character_),
    stringsAsFactors = FALSE)

fxGtoPdb <- function() data.frame(
    QueryIDs = c("imatinib", "imatinib"),
    ligandName = c("imatinib", "imatinib"),
    target_gene = c("ABL1", "DDR1"),
    type = c("Inhibitor", "Inhibitor"),
    action = c("Inhibition", "Inhibition"),
    refIds = c("12,34", NA_character_),
    stringsAsFactors = FALSE)

fxDgidb <- function() data.frame(
    QueryIDs = c("imatinib", "imatinib"),
    drug_name = c("IMATINIB", "IMATINIB"),
    drug_aliases = rep("CHEMBL:CHEMBL941; DRUGBANK:DB00619", 2),
    gene_name = c("KIT", "MAPK10"),
    interaction_types = c("inhibitor", ""),   ## DGIdb often has none
    sources = c("TEND; FDA", "DTC"),
    stringsAsFactors = FALSE)

## --- schema and basic behaviour ---------------------------------------

test_that("assembleMoaTable returns the documented columns, even when empty", {
    empty <- assembleMoaTable(list())
    expect_s3_class(empty, "data.frame")
    expect_identical(nrow(empty), 0L)
    expect_identical(names(empty), .dtiMoaCols)

    out <- assembleMoaTable(list(chembl = fxChembl()))
    expect_identical(names(out), .dtiMoaCols)
    expect_identical(nrow(out), 2L)
})

test_that("assembleMoaTable rejects sources it has no MOA mapping for", {
    expect_error(assembleMoaTable(list(pubchem = data.frame(x = 1))),
                 "no MOA mapping")
})

test_that("a source present but matching zero rows contributes nothing", {
    zero <- fxChembl()[0, , drop = FALSE]
    out <- assembleMoaTable(list(chembl = zero, opentargets = fxOpenTargets()))
    expect_identical(unique(out$source), "opentargets")

    ## and all-empty is the empty table, not an error
    expect_identical(nrow(assembleMoaTable(list(chembl = zero))), 0L)
})

test_that("rows carrying no mechanism information at all are dropped", {
    ## TTD's MCL1 row has MOA = NA; DGIdb's MAPK10 row has "" -> both go.
    ttd <- assembleMoaTable(list(ttd = fxTtd()))
    expect_identical(nrow(ttd), 1L)
    expect_identical(ttd$target_symbol, "BCR-ABL1")

    dgidb <- assembleMoaTable(list(dgidb = fxDgidb()))
    expect_identical(nrow(dgidb), 1L)
    expect_identical(dgidb$target_symbol, "KIT")
})

## --- grain: the invariant this whole layer exists to protect ----------

test_that("Broad's drug-level MOA never carries a target attribution", {
    out <- assembleMoaTable(list(broad = fxBroad()))
    expect_true(all(out$moa_scope == "drug"))
    ## despite fxBroad() having a populated target_gene column
    expect_true(all(is.na(out$target_symbol)))
    expect_true(all(is.na(out$target_uniprot)))
})

test_that("Broad's ' | '-joined MOA field explodes to one row per mechanism", {
    out <- assembleMoaTable(list(broad = fxBroad()))
    ## 3 identical target rows x 2 mechanisms, deduped to the 2 distinct
    ## drug-level statements - the repetition across target rows carries
    ## no extra information.
    expect_identical(sort(out$moa_text),
                     c("Bcr-Abl kinase inhibitor", "KIT inhibitor"))
    expect_identical(nrow(out), 2L)
})

test_that("pair-level sources keep their target attribution", {
    out <- assembleMoaTable(list(chembl = fxChembl(), opentargets = fxOpenTargets()))
    expect_true(all(out$moa_scope == "drug_target"))
    expect_identical(sort(out$target_uniprot[out$source == "chembl"]),
                     c("P00519", "P10721"))
    expect_identical(sort(out$target_symbol[out$source == "opentargets"]),
                     c("ABL1", "BCR"))
})

test_that("one Open Targets mechanism can span several targets", {
    ## The discriminator for finding (1): a single MOA string legitimately
    ## maps to >1 distinct target, so MOA text is not a per-target key.
    out <- assembleMoaTable(list(opentargets = fxOpenTargets()))
    byMoa <- split(out$target_symbol, out$moa_text)
    expect_true(any(vapply(byMoa, function(t) length(unique(t)) > 1L, logical(1))))
})

test_that("the intended key does not multiply rows", {
    out <- assembleMoaTable(list(chembl = fxChembl(), opentargets = fxOpenTargets(),
                                 broad = fxBroad(), ttd = fxTtd(),
                                 gtopdb = fxGtoPdb(), dgidb = fxDgidb()))
    ## Both target columns belong in the key: a source may name its target
    ## by accession only (ChEMBL) or by symbol only (most others), so
    ## keying on target_symbol alone collides two genuinely distinct
    ## ChEMBL targets that share one mechanism string.
    key <- paste(out$drug_key, out$source, out$moa_text, out$target_symbol,
                 out$target_uniprot, sep = "\r")
    expect_identical(anyDuplicated(key), 0L)
})

test_that("one mechanism over two accessions stays two rows, not a collision", {
    ## ChEMBL reports "Bcr/Abl fusion protein inhibitor" against both
    ## P00519 (ABL1) and P11274 (BCR) - the same fan-out Open Targets
    ## shows by symbol. These must survive as distinct rows.
    fx <- data.frame(
        QueryIDs = c("CHEMBL941", "CHEMBL941"),
        chembl_id = c("CHEMBL941", "CHEMBL941"),
        Drug_Name = c("IMATINIB", "IMATINIB"),
        MOA = rep("Bcr/Abl fusion protein inhibitor", 2),
        Action_Type = c("INHIBITOR", "INHIBITOR"),
        UniProt_ID = c("P00519", "P11274"),
        stringsAsFactors = FALSE)
    out <- assembleMoaTable(list(chembl = fx))
    expect_identical(nrow(out), 2L)
    expect_identical(sort(out$target_uniprot), c("P00519", "P11274"))
})

test_that("resolveGeneSymbol is off by default and never invents a symbol", {
    ## Default path must not touch the network, so a bad accession is
    ## simply left alone rather than looked up.
    out <- assembleMoaTable(list(chembl = fxChembl()))
    expect_true(all(is.na(out$target_symbol)))
    expect_identical(sort(out$target_uniprot), c("P00519", "P10721"))
})

## --- kind, normalization, keys ----------------------------------------

test_that("moa_kind separates free-text prose from controlled vocabulary", {
    out <- assembleMoaTable(list(chembl = fxChembl(), ttd = fxTtd()))
    expect_identical(unique(out$moa_kind[out$source == "chembl"]), "description")
    expect_identical(unique(out$moa_kind[out$source == "ttd"]), "action_vocabulary")
})

test_that("action terms normalize across sources' differing case and wording", {
    out <- assembleMoaTable(list(chembl = fxChembl(), ttd = fxTtd(),
                                 gtopdb = fxGtoPdb(), dgidb = fxDgidb()))
    ## ChEMBL "INHIBITOR", TTD "Inhibitor", GtoPdb "Inhibitor", DGIdb
    ## "inhibitor" all land on one term while raw values stay verbatim.
    expect_identical(unique(out$action_type), "inhibitor")
    expect_true(all(c("INHIBITOR", "Inhibitor", "inhibitor") %in%
                    out$action_type_raw))
})

test_that("action normalization keeps near-miss terms distinct", {
    ## "antagonist" contains "agonist"; the two inverse/partial forms must
    ## not collapse into plain "agonist" either.
    expect_identical(
        .dtiNormalizeActionType(c("ANTAGONIST", "Agonist", "inverse agonist",
                                  "PARTIAL AGONIST", "Inhibition",
                                  "Channel blocker", "Allosteric modulator")),
        c("antagonist", "agonist", "inverse agonist", "partial agonist",
          "inhibitor", "blocker", "modulator"))
})

test_that("unclassifiable action terms are NA rather than guessed", {
    expect_true(all(is.na(.dtiNormalizeActionType(c("CAR-T-Cell-Therapy", "", NA)))))
})

test_that("drug_key prefers a ChEMBL ID, then InChIKey, then folded name", {
    out <- assembleMoaTable(list(chembl = fxChembl(), opentargets = fxOpenTargets(),
                                 broad = fxBroad(), ttd = fxTtd(),
                                 dgidb = fxDgidb()))
    keyOf <- function(s) unique(out$drug_key[out$source == s])
    expect_identical(keyOf("chembl"), "CHEMBL941")
    expect_identical(keyOf("opentargets"), "CHEMBL941")
    expect_identical(keyOf("dgidb"), "CHEMBL941")   ## mined from drug_aliases
    expect_identical(keyOf("broad"), "KTUFNOKKBVMGRW-UHFFFAOYSA-N")
    expect_identical(keyOf("ttd"), "imatinib")      ## name fallback, folded
})

test_that("query_id uses the caller's original token when 'resolved' is present", {
    res <- list(chembl = fxChembl())
    ## queryDrugTargets() records original -> per-source-resolved
    attr(res, "resolved") <- list(chembl = c(imatinib = "CHEMBL941"))
    expect_identical(unique(assembleMoaTable(res)$query_id), "imatinib")
    ## and falls back to the source's own QueryIDs otherwise
    expect_identical(unique(assembleMoaTable(list(chembl = fxChembl()))$query_id),
                     "CHEMBL941")
})

test_that("evidence_refs carries per-source provenance where the source has it", {
    out <- assembleMoaTable(list(gtopdb = fxGtoPdb(), dgidb = fxDgidb()))
    expect_identical(out$evidence_refs[out$source == "gtopdb" &
                                       out$target_symbol == "ABL1"], "12,34")
    expect_identical(out$evidence_refs[out$source == "dgidb"], "TEND; FDA")
})

## --- reshapers --------------------------------------------------------

test_that("moaWide keeps every distinct MOA string (round trip)", {
    long <- assembleMoaTable(list(chembl = fxChembl(), opentargets = fxOpenTargets(),
                                  ttd = fxTtd()))
    wide <- moaWide(long)
    expect_identical(nrow(wide), length(unique(paste(long$query_id, long$drug_key))))
    flat <- unlist(c(wide$moa_description, wide$moa_action_vocabulary))
    expect_setequal(flat, unique(long$moa_text[!is.na(long$moa_text)]))
})

test_that("moaWide separates the two MOA kinds and counts sources", {
    wide <- moaWide(assembleMoaTable(list(chembl = fxChembl(), ttd = fxTtd())))
    ## fxChembl and fxTtd key to different drug_keys (CHEMBL941 vs the name
    ## fallback), so one row each - both carry only their own kind.
    chemblRow <- wide[wide$drug_key == "CHEMBL941", ]
    expect_length(chemblRow$moa_action_vocabulary[[1L]], 0L)
    expect_true(length(chemblRow$moa_description[[1L]]) > 0L)
    expect_identical(chemblRow$n_sources, 1L)
})

test_that("moaWide on an empty table returns the empty wide shape", {
    wide <- moaWide(assembleMoaTable(list()))
    expect_identical(nrow(wide), 0L)
    expect_true(all(c("query_id", "drug_key", "n_sources", "sources",
                      "moa_description", "moa_action_vocabulary",
                      "action_type") %in% names(wide)))
})

test_that("moaByTarget drops drug-scope rows and keeps the rest", {
    long <- assembleMoaTable(list(chembl = fxChembl(), broad = fxBroad()))
    expect_true(any(long$moa_scope == "drug"))
    expect_message(moaByTarget(long), "drug-level MOA row")
    expect_silent(byTgt <- moaByTarget(long, verbose = FALSE))
    expect_false(any(byTgt$moa_scope == "drug"))
    expect_identical(unique(byTgt$source), "chembl")
})

test_that("the reshapers reject anything that is not a MOA table", {
    expect_error(moaWide(data.frame(a = 1)), "missing column")
    expect_error(moaByTarget("nope"), "Expected a data.frame")
})

## --- source capability reference --------------------------------------

test_that("listMoaSources documents every source, including the excluded one", {
    ms <- listMoaSources()
    expect_setequal(ms$source, c(names(.dtiMoaSourceSpec), "pubchem"))
    expect_true(is.na(ms$moa_column[ms$source == "pubchem"]))
    ## the spec table and the reference table must not drift apart
    for (src in names(.dtiMoaSourceSpec)) {
        expect_identical(ms$moa_column[ms$source == src],
                         .dtiMoaSourceSpec[[src]]$moa)
        expect_identical(ms$moa_scope[ms$source == src],
                         .dtiMoaSourceSpec[[src]]$scope)
    }
})

## --- live known-answer tests ------------------------------------------

test_that("live: ChEMBL returns imatinib's mechanisms despite the salt filing", {
    skip_if_offline_dti()
    ## Regression guard for the parent/salt fix: ChEMBL files imatinib's
    ## mechanism records against the mesylate salt CHEMBL1642, so querying
    ## molecule_chembl_id alone returned zero rows.
    res <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL941"))
    expect_gt(nrow(res), 0L)
    expect_identical(unique(res$chembl_id), "CHEMBL941")
    expect_identical(unique(res$Drug_Name), "IMATINIB")
    expect_true("Tyrosine-protein kinase ABL inhibitor" %in% res$MOA)

    moa <- assembleMoaTable(list(chembl = res))
    expect_true("Tyrosine-protein kinase ABL inhibitor" %in% moa$moa_text)
    expect_true(all(moa$moa_scope == "drug_target"))
})

test_that("live: resolveGeneSymbol aligns ChEMBL's accessions with other sources", {
    skip_if_offline_dti()
    ## Without it, ChEMBL rows carry only target_uniprot and cannot be
    ## grouped per target alongside symbol-keyed sources.
    res <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL941"))
    plain <- assembleMoaTable(list(chembl = res))
    expect_true(all(is.na(plain$target_symbol)))

    filled <- assembleMoaTable(list(chembl = res), resolveGeneSymbol = TRUE)
    expect_true(any(!is.na(filled$target_symbol)))
    expect_identical(filled$target_symbol[filled$target_uniprot == "P00519"][1L],
                     "ABL1")
    ## fills only the gap; accessions and every other column are untouched
    expect_identical(filled$target_uniprot, plain$target_uniprot)
    expect_identical(filled$moa_text, plain$moa_text)
})

test_that("live: dasatinib's ChEMBL MOA survives assembly verbatim", {
    skip_if_offline_dti()
    res <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL1421"))
    moa <- assembleMoaTable(list(chembl = res))
    expect_true("Tyrosine-protein kinase ABL inhibitor" %in% moa$moa_text)
    expect_identical(unique(moa$moa_kind), "description")
})

test_that("live: Open Targets fans one mechanism across several targets", {
    skip_if_offline_dti()
    res <- getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
                                          ids = "imatinib"))
    moa <- assembleMoaTable(list(opentargets = res))
    expect_true("Bcr/Abl fusion protein inhibitor" %in% moa$moa_text)
    spanning <- split(moa$target_symbol, moa$moa_text)
    expect_true(any(vapply(spanning, function(t) length(unique(t)) > 1L, logical(1))))
})

test_that("live: TTD contributes controlled-vocabulary MOA, not prose", {
    skip_if_offline_dti()
    dbPath <- buildTtdDb(rerun = FALSE)
    res <- ttdTargetAnnot(list(molType = "cmp", idType = "name", ids = "Imatinib"),
                          dbPath)
    moa <- assembleMoaTable(list(ttd = res))
    expect_identical(unique(moa$moa_kind), "action_vocabulary")
    expect_true("Inhibitor" %in% moa$moa_text)
    expect_identical(unique(moa$action_type), "inhibitor")
})
