## Tests for the cross-source MOA layer: assembleMoaTable() (drug -> MOA)
## and assembleMoaTargets() (MOA -> target link).
##
## Both are pure transformers over queryDrugTargets()'s output, so most
## coverage here is network-free and runs against synthetic per-source
## frames - deliberate, since it exercises the concept boundaries that
## matter (MOA is drug-level; only three sources carry MOA at all)
## without needing several live APIs up at once.
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
## Column names mirror what each accessor really returns, trimmed to what
## the extractor reads. fxChembl() encodes the central case: ONE mechanism
## of imatinib recorded against TWO targets.

fxChembl <- function() data.frame(
    QueryIDs = c("CHEMBL941", "CHEMBL941", "CHEMBL941"),
    chembl_id = "CHEMBL941",
    Drug_Name = "IMATINIB",
    MOA = c("Bcr/Abl fusion protein inhibitor",
            "Bcr/Abl fusion protein inhibitor",
            "Stem cell growth factor receptor inhibitor"),
    Action_Type = "INHIBITOR",
    UniProt_ID = c("P00519", "P11274", "P10721"),
    stringsAsFactors = FALSE)

fxOpenTargets <- function() data.frame(
    QueryIDs = c("CHEMBL941", "CHEMBL941"),
    chembl_id = "CHEMBL941",
    drug_name = "IMATINIB",
    mechanism_of_action = "Bcr/Abl fusion protein inhibitor",
    action_type = "INHIBITOR",
    approved_symbol = c("ABL1", "BCR"),
    stringsAsFactors = FALSE)

## Broad repeats one " | "-joined drug-level string on every target row.
fxBroad <- function() data.frame(
    QueryIDs = "imatinib",
    pert_iname = "imatinib",
    InChIKey = "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
    moa = "Bcr-Abl kinase inhibitor | KIT inhibitor",
    target_gene = c("BCR", "ABL1", "CSF1R"),
    stringsAsFactors = FALSE)

## --- source selection: the concept boundary ---------------------------

test_that("only the three MOA-carrying sources are accepted", {
    for (src in c("ttd", "gtopdb", "dgidb", "pubchem")) {
        res <- stats::setNames(list(data.frame(x = 1)), src)
        expect_error(assembleMoaTable(res), "does not use source")
        expect_error(assembleMoaTargets(res), "does not use source")
    }
})

test_that("the rejection points at where that data actually lives", {
    expect_error(assembleMoaTable(list(ttd = data.frame(x = 1))),
                 "combineDrugTargets")
})

test_that("listMoaSources documents all seven sources and the exclusions", {
    ms <- listMoaSources()
    expect_setequal(ms$source, c("chembl", "opentargets", "broad", "ttd",
                                 "gtopdb", "dgidb", "pubchem"))
    expect_setequal(ms$source[ms$used_for_moa], names(.dtiMoaSourceSpec))
    ## the three excluded action-type sources are labelled as such
    expect_identical(sort(ms$source[ms$kind == "action type"]),
                     c("dgidb", "gtopdb", "ttd"))
    ## the cardinality evidence that separates the two concepts
    expect_lt(ms$n_distinct[ms$source == "ttd"],
              ms$n_distinct[ms$source == "chembl"])
})

## --- MOA is drug-level: the central property --------------------------

test_that("one mechanism against two targets is ONE MOA row", {
    ## The whole reason the two relations are separate: ChEMBL records
    ## "Bcr/Abl fusion protein inhibitor" twice (P00519, P11274), but the
    ## drug has one such mechanism, not two.
    moa <- assembleMoaTable(list(chembl = fxChembl()))
    expect_identical(nrow(moa), 2L)
    expect_setequal(moa$moa_text,
                    c("Bcr/Abl fusion protein inhibitor",
                      "Stem cell growth factor receptor inhibitor"))
})

test_that("the same mechanism keeps both of its targets in the link table", {
    lnk <- assembleMoaTargets(list(chembl = fxChembl()))
    expect_identical(nrow(lnk), 3L)
    bcrabl <- lnk[lnk$moa_text == "Bcr/Abl fusion protein inhibitor", ]
    expect_setequal(bcrabl$target_uniprot, c("P00519", "P11274"))
})

test_that("MOA table carries no target columns at all", {
    ## Targets are a separate relation; leaking them here would re-create
    ## the per-edge framing the two-table model exists to avoid.
    expect_false(any(c("target_symbol", "target_uniprot") %in%
                     names(assembleMoaTable(list(chembl = fxChembl())))))
})

test_that("a drug can carry several distinct MOA terms", {
    moa <- assembleMoaTable(list(chembl = fxChembl()))
    expect_identical(length(unique(moa$moa_text)), 2L)
    expect_identical(unique(moa$drug_key), "CHEMBL941")
})

## --- Broad: states MOA, links no targets ------------------------------

test_that("Broad contributes MOA terms but no link rows", {
    both <- list(chembl = fxChembl(), broad = fxBroad())
    moa <- assembleMoaTable(both)
    expect_true("broad" %in% moa$source)
    expect_setequal(moa$moa_text[moa$source == "broad"],
                    c("Bcr-Abl kinase inhibitor", "KIT inhibitor"))
    ## despite fxBroad() having a populated target_gene column
    expect_false("broad" %in% assembleMoaTargets(both)$source)
})

test_that("Broad's ' | '-joined field explodes to one row per MOA term", {
    moa <- assembleMoaTable(list(broad = fxBroad()))
    ## 3 identical target rows x 2 terms -> the 2 distinct drug-level terms
    expect_identical(nrow(moa), 2L)
})

## --- schema, dedup, empties -------------------------------------------

test_that("both assemblers return their documented columns, even when empty", {
    expect_identical(names(assembleMoaTable(list())), .dtiMoaCols)
    expect_identical(names(assembleMoaTargets(list())), .dtiMoaTargetCols)
    expect_identical(nrow(assembleMoaTable(list())), 0L)
    expect_identical(names(assembleMoaTable(list(chembl = fxChembl()))),
                     .dtiMoaCols)
    expect_identical(names(assembleMoaTargets(list(chembl = fxChembl()))),
                     .dtiMoaTargetCols)
})

test_that("a source present but matching zero rows contributes nothing", {
    zero <- fxChembl()[0, , drop = FALSE]
    out <- assembleMoaTable(list(chembl = zero, broad = fxBroad()))
    expect_identical(unique(out$source), "broad")
    expect_identical(nrow(assembleMoaTable(list(chembl = zero))), 0L)
    expect_identical(nrow(assembleMoaTargets(list(chembl = zero))), 0L)
})

test_that("rows with no MOA term are dropped", {
    fx <- fxChembl(); fx$MOA <- c("Bcr/Abl fusion protein inhibitor", NA, "")
    moa <- assembleMoaTable(list(chembl = fx))
    expect_identical(nrow(moa), 1L)
})

test_that("neither table multiplies rows", {
    both <- list(chembl = fxChembl(), opentargets = fxOpenTargets(),
                 broad = fxBroad())
    moa <- assembleMoaTable(both)
    expect_identical(anyDuplicated(
        paste(moa$drug_key, moa$source, moa$moa_text, sep = "\r")), 0L)
    lnk <- assembleMoaTargets(both)
    expect_identical(anyDuplicated(
        paste(lnk$drug_key, lnk$source, lnk$moa_text, lnk$target_symbol,
              lnk$target_uniprot, sep = "\r")), 0L)
})

test_that("sources are never merged with each other", {
    ## ChEMBL and Open Targets state the same mechanism; merging them
    ## would be a cross-source identity claim this package does not make.
    moa <- assembleMoaTable(list(chembl = fxChembl(),
                                 opentargets = fxOpenTargets()))
    shared <- moa[moa$moa_text == "Bcr/Abl fusion protein inhibitor", ]
    expect_setequal(shared$source, c("chembl", "opentargets"))
})

## --- action types and keys --------------------------------------------

test_that("action terms normalize while the raw term is kept", {
    moa <- assembleMoaTable(list(chembl = fxChembl()))
    expect_identical(unique(moa$action_type), "inhibitor")
    expect_identical(unique(moa$action_type_raw), "INHIBITOR")
})

test_that("action normalization keeps near-miss terms distinct", {
    ## "antagonist" contains "agonist"; inverse/partial must not collapse.
    expect_identical(
        .dtiNormalizeActionType(c("ANTAGONIST", "Agonist", "inverse agonist",
                                  "PARTIAL AGONIST", "Inhibition",
                                  "Channel blocker", "Allosteric modulator")),
        c("antagonist", "agonist", "inverse agonist", "partial agonist",
          "inhibitor", "blocker", "modulator"))
})

test_that("unclassifiable or absent action terms are NA, not guessed", {
    expect_true(all(is.na(.dtiNormalizeActionType(c("Radiotherapy agent", "", NA)))))
    ## Broad states no action type at all
    moa <- assembleMoaTable(list(broad = fxBroad()))
    expect_true(all(is.na(moa$action_type_raw)))
})

test_that("drug_key prefers a ChEMBL ID, then InChIKey, then folded name", {
    moa <- assembleMoaTable(list(chembl = fxChembl(),
                                 opentargets = fxOpenTargets(),
                                 broad = fxBroad()))
    keyOf <- function(s) unique(moa$drug_key[moa$source == s])
    expect_identical(keyOf("chembl"), "CHEMBL941")
    expect_identical(keyOf("opentargets"), "CHEMBL941")
    expect_identical(keyOf("broad"), "KTUFNOKKBVMGRW-UHFFFAOYSA-N")

    ## name fallback, case-folded, when no structured key is available
    noKey <- fxBroad(); noKey$InChIKey <- NA_character_
    expect_identical(unique(assembleMoaTable(list(broad = noKey))$drug_key),
                     "imatinib")
})

test_that("query_id uses the caller's original token when 'resolved' is present", {
    res <- list(chembl = fxChembl())
    attr(res, "resolved") <- list(chembl = c(imatinib = "CHEMBL941"))
    expect_identical(unique(assembleMoaTable(res)$query_id), "imatinib")
    expect_identical(unique(assembleMoaTable(list(chembl = fxChembl()))$query_id),
                     "CHEMBL941")
})

test_that("mechanism_comment rides along when ChEMBL was fetched with fields='all'", {
    fx <- fxChembl()
    fx$mechanism.mechanism_comment <- c("Role in regulating gastric secretion",
                                        NA, NA)
    lnk <- assembleMoaTargets(list(chembl = fx))
    expect_identical(lnk$mechanism_comment[lnk$target_uniprot == "P00519"],
                     "Role in regulating gastric secretion")
    ## absent under the default column set, without erroring
    expect_true(all(is.na(assembleMoaTargets(list(chembl = fxChembl()))$mechanism_comment)))
})

## --- reshapers and validation -----------------------------------------

test_that("moaWide keeps every distinct MOA term (round trip)", {
    long <- assembleMoaTable(list(chembl = fxChembl(), broad = fxBroad()))
    wide <- moaWide(long)
    expect_identical(nrow(wide),
                     length(unique(paste(long$query_id, long$drug_key))))
    expect_setequal(unlist(wide$moa_text), unique(long$moa_text))
    expect_identical(sum(wide$n_moa), length(unique(paste(long$drug_key,
                                                          long$moa_text))))
})

test_that("moaWide on an empty table returns the empty wide shape", {
    wide <- moaWide(assembleMoaTable(list()))
    expect_identical(nrow(wide), 0L)
    expect_true(all(c("query_id", "drug_key", "n_sources", "n_moa", "sources",
                      "moa_text", "action_type") %in% names(wide)))
})

test_that("moaWide rejects anything that is not a MOA table", {
    expect_error(moaWide(data.frame(a = 1)), "missing column")
    expect_error(moaWide("nope"), "Expected a data.frame")
    ## the link table is not a MOA table
    expect_error(moaWide(assembleMoaTargets(list(chembl = fxChembl()))),
                 "missing column")
})

## --- live known-answer tests ------------------------------------------

test_that("live: ChEMBL's imatinib mechanisms survive the salt filing", {
    skip_if_offline_dti()
    ## Regression guard: ChEMBL files imatinib's mechanism records against
    ## the mesylate salt CHEMBL1642, so molecule_chembl_id alone gave none.
    res <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL941"))
    expect_gt(nrow(res), 0L)
    expect_identical(unique(res$Drug_Name), "IMATINIB")

    moa <- assembleMoaTable(list(chembl = res))
    expect_true("Tyrosine-protein kinase ABL inhibitor" %in% moa$moa_text)
    ## MOA rows collapse the target fan-out; link rows do not
    expect_lt(nrow(moa), nrow(assembleMoaTargets(list(chembl = res))))
})

test_that("live: a real MOA term is reused across different drugs", {
    skip_if_offline_dti()
    ## The property that makes MOA a drug-level descriptor rather than a
    ## per-edge fact: the same term describes unrelated drugs.
    a <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                  ids = "CHEMBL25"))          # aspirin
    b <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                  ids = "CHEMBL521"))         # ibuprofen
    moa <- assembleMoaTable(list(chembl = rbind(a, b)))
    shared <- intersect(moa$moa_text[moa$drug_key == "CHEMBL25"],
                        moa$moa_text[moa$drug_key == "CHEMBL521"])
    expect_true(length(shared) > 0L)
})

test_that("live: Open Targets states the same mechanism over several targets", {
    skip_if_offline_dti()
    res <- getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
                                          ids = "imatinib"))
    moa <- assembleMoaTable(list(opentargets = res))
    lnk <- assembleMoaTargets(list(opentargets = res))
    expect_true("Bcr/Abl fusion protein inhibitor" %in% moa$moa_text)
    spanning <- lnk$target_symbol[lnk$moa_text == "Bcr/Abl fusion protein inhibitor"]
    expect_gt(length(unique(spanning)), 1L)
})

test_that("live: resolveGeneSymbol aligns ChEMBL accessions with symbols", {
    skip_if_offline_dti()
    res <- getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
                                    ids = "CHEMBL941"))
    plain <- assembleMoaTargets(list(chembl = res))
    expect_true(all(is.na(plain$target_symbol)))

    filled <- assembleMoaTargets(list(chembl = res), resolveGeneSymbol = TRUE)
    expect_identical(filled$target_symbol[filled$target_uniprot == "P00519"][1L],
                     "ABL1")
    expect_identical(filled$target_uniprot, plain$target_uniprot)
    expect_identical(filled$moa_text, plain$moa_text)
})


## --- buildMoaMasterTable() -------------------------------------------
## The Broad half reads a local SQLite, so it is exercised against a
## synthetic one and needs no network. Only the ChEMBL sweep is guarded.

.moaSyntheticBroadDb <- function() {
    path <- tempfile(fileext = ".db")
    con <- RSQLite::dbConnect(RSQLite::SQLite(), path)
    on.exit(RSQLite::dbDisconnect(con))
    RSQLite::dbWriteTable(con, "broad_interactions", data.frame(
        pert_iname = c("drugA", "drugA", "drugB", "drugC", "drugD"),
        ## drugA is annotated against two genes; drugB packs two
        ## mechanisms into one string; drugC has a mechanism but no
        ## target at all; drugD has a target but no mechanism.
        target_gene = c("FGFR1", "KLB", "", NA_character_, "ABL1"),
        moa = c("kinase inhibitor", "kinase inhibitor",
                "dopamine receptor antagonist | serotonin receptor antagonist",
                "antitumor agent", ""),
        stringsAsFactors = FALSE))
    path
}

test_that("buildMoaMasterTable enumerates Broad drugs, one row per mechanism", {
    out <- buildMoaMasterTable(sources = "broad",
                               brhDbPath = .moaSyntheticBroadDb(), verbose = FALSE)
    expect_identical(names(out), c("drug_id", "drug_name", "moa", "action",
                                   "source", "has_target"))
    expect_true(all(out$source == "broad"))
    ## drugD has no mechanism, so it is not a row here at all.
    expect_setequal(unique(out$drug_id), c("drugA", "drugB", "drugC"))
    ## drugA's two target rows collapse to one mechanism row.
    expect_identical(sum(out$drug_id == "drugA"), 1L)
    ## drugB's packed string becomes two rows, one mechanism each.
    expect_setequal(out$moa[out$drug_id == "drugB"],
                    c("dopamine receptor antagonist", "serotonin receptor antagonist"))
    expect_false(any(grepl("|", out$moa, fixed = TRUE)))
    ## Broad records no action type of its own.
    expect_true(all(is.na(out$action)))
})

test_that("buildMoaMasterTable reports whether a target was named", {
    out <- buildMoaMasterTable(sources = "broad",
                               brhDbPath = .moaSyntheticBroadDb(), verbose = FALSE)
    expect_true(out$has_target[out$drug_id == "drugA"])
    ## Empty string and NA both count as no target named.
    expect_false(any(out$has_target[out$drug_id == "drugB"]))
    expect_false(out$has_target[out$drug_id == "drugC"])
})

test_that("buildMoaMasterTable holds one row per source, drug and mechanism", {
    out <- buildMoaMasterTable(sources = "broad",
                               brhDbPath = .moaSyntheticBroadDb(), verbose = FALSE)
    expect_false(any(duplicated(out[, c("source", "drug_id", "moa")])))
    expect_false(any(is.na(out$moa) | !nzchar(out$moa)))
})

test_that("buildMoaMasterTable validates its arguments", {
    expect_error(buildMoaMasterTable(sources = "Broad"), 'Did you mean "broad"')
    expect_error(buildMoaMasterTable(sources = "opentargets"),
                 "does not recognise: opentargets")
    expect_error(buildMoaMasterTable(sources = "broad"), "brhDbPath")
})

test_that("buildMoaMasterTable sweeps ChEMBL's mechanism collection", {
    skip_if_offline_dti()
    out <- buildMoaMasterTable(sources = "chembl", resolveDrugNames = FALSE,
                               verbose = FALSE)
    expect_gt(nrow(out), 5000L)
    expect_true(all(out$source == "chembl"))
    expect_true(all(grepl("^CHEMBL", out$drug_id)))
    ## "Unknown" is a placeholder, not a mechanism.
    expect_false(any(tolower(out$moa) == "unknown"))
    expect_true(any(tolower(buildMoaMasterTable(
        sources = "chembl", resolveDrugNames = FALSE, includeUnknown = TRUE,
        verbose = FALSE)$moa) == "unknown"))
    ## The population this table exists to reach: annotated, untargeted.
    expect_true(any(!out$has_target))
    ## Known answer: methazolamide is a carbonic anhydrase inhibitor.
    expect_true(any(out$drug_id == "CHEMBL19" &
                    grepl("Carbonic anhydrase", out$moa)))
})
