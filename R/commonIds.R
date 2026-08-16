## ---------------------------------------------------------------------
##  commonIds.R - shared join keys across the per-source tables
##
##  Each source returns its own identifier vocabulary: ChEMBL is keyed on
##  UniProt accessions, Open Targets on Ensembl gene IDs, DGIdb/Broad/GtoPdb
##  on gene symbols, TTD on both an accession and a symbol. Nothing joins
##  without translation, and row-appending them produces a table whose key
##  columns are differently named and differently populated per block.
##
##  This file adds four canonical columns to those tables - hgnc_id,
##  gene_symbol, target_uniprot, compound_chembl_id - so they can be merged
##  or stacked. Two design rules matter:
##
##  1. Lookups, never joins. Every column is filled with match() against a
##     many-to-one lookup vector, so the row count cannot change. A merge()
##     on HGNC would silently multiply rows for the 42 genes carrying more
##     than one UniProt accession, which is exactly the grain corruption
##     this layer must not introduce.
##  2. Ambiguity yields NA, not a guess. 116 UniProt accessions in HGNC map
##     to more than one gene, and 42 genes map to more than one accession.
##     In both directions the ambiguous cases are left NA and hgnc_id stays
##     the reliable key, rather than picking an arbitrary winner that would
##     look authoritative and be wrong.
##
##  Compound-side coverage is currently ChEMBL IDs only, which ChEMBL, Open
##  Targets and part of DGIdb already carry natively. TTD, Broad and GtoPdb
##  key their compounds on PubChem CIDs, InChIKeys and IUPHAR ligand IDs
##  respectively; mapping those to ChEMBL needs UniChem and is a separate,
##  more expensive phase.
## ---------------------------------------------------------------------

#' Per-source columns the canonical keys are derived from
#'
#' \code{NA} means the source has no column of that kind at all. The
#' compound column is the source's own drug identifier; only the entries
#' flagged in \code{.dtiChemblIdNative} already hold a ChEMBL ID.
#' @keywords internal
#' @noRd
.dtiCommonIdSpec <- list(
    chembl      = list(uniprot = "UniProt_ID", symbol = NA_character_,
                       ensembl = NA_character_, compound = "chembl_id"),
    dgidb       = list(uniprot = NA_character_, symbol = "gene_name",
                       ensembl = NA_character_, compound = "drug_concept_id"),
    opentargets = list(uniprot = NA_character_, symbol = "approved_symbol",
                       ensembl = "ensembl_id", compound = "drug_id"),
    ttd         = list(uniprot = "Uniprot_acc", symbol = "GeneName",
                       ensembl = NA_character_, compound = NA_character_),
    broad       = list(uniprot = NA_character_, symbol = "target_gene",
                       ensembl = NA_character_, compound = NA_character_),
    gtopdb      = list(uniprot = NA_character_, symbol = "target_gene",
                       ensembl = NA_character_, compound = NA_character_)
)

#' Canonical columns this layer adds, in output order
#' @keywords internal
#' @noRd
.dtiCommonIdCols <- c("hgnc_id", "gene_symbol", "target_uniprot",
                      "compound_chembl_id")

#' Build the many-to-one lookup vectors used to fill the canonical columns
#'
#' Each returned vector is named by the foreign identifier and holds the
#' HGNC id (or accession) it resolves to. Identifiers resolving to more
#' than one gene are dropped rather than resolved arbitrarily.
#' @keywords internal
#' @noRd
.dtiHgncLookups <- function(hgncTable) {
    uniqueOnly <- function(keys, values) {
        keep <- !is.na(keys) & nzchar(keys)
        keys <- keys[keep]; values <- values[keep]
        dup <- unique(keys[duplicated(keys)])
        keep <- !keys %in% dup
        stats::setNames(values[keep], keys[keep])
    }
    acc <- .hgncExplode(hgncTable, "uniprot_ids")
    single <- lengths(hgncTable$uniprot_ids) == 1L

    list(
        bySymbol  = stats::setNames(hgncTable$hgnc_id, hgncTable$symbol),
        byAcc     = uniqueOnly(acc$id, acc$hgnc_id),
        byEnsembl = uniqueOnly(hgncTable$ensembl_gene_id, hgncTable$hgnc_id),
        ## hgnc_id -> its single accession; genes with several are absent,
        ## so they resolve to NA rather than to whichever came first.
        accOf     = uniqueOnly(hgncTable$hgnc_id[single],
                               vapply(hgncTable$uniprot_ids[single], `[`, character(1), 1L)),
        symbolOf  = stats::setNames(hgncTable$symbol, hgncTable$hgnc_id)
    )
}

#' Pull a ChEMBL molecule ID out of a source's own compound identifier
#'
#' ChEMBL and Open Targets hold a bare ChEMBL ID; DGIdb holds a prefixed
#' concept id that is a ChEMBL ID only when the prefix says so
#' (\code{"chembl:CHEMBL52885"}), and something else entirely otherwise
#' (\code{"ncit:C104267"}, \code{"rxcui:357977"}).
#' @keywords internal
#' @noRd
.dtiAsChemblId <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    out <- ifelse(grepl("^chembl:", x, ignore.case = TRUE), sub("^[Cc][Hh][Ee][Mm][Bb][Ll]:", "", x),
                  ifelse(grepl("^CHEMBL[0-9]+$", x), x, NA_character_))
    out
}

#' Add shared gene and compound identifiers to per-source result tables
#'
#' The source functions each return their own identifier vocabulary -
#' ChEMBL uses UniProt accessions, Open Targets Ensembl gene IDs, DGIdb,
#' the Broad Repurposing Hub and GtoPdb gene symbols - which leaves no
#' column to merge or stack them on. This function adds four columns to
#' every table so they share join keys: \code{hgnc_id}, \code{gene_symbol}
#' (the current approved symbol), \code{target_uniprot} and
#' \code{compound_chembl_id}.
#'
#' Gene identifiers are resolved against the HGNC gene table, which maps
#' symbols, UniProt accessions and Ensembl gene IDs to one another
#' offline, so no network access is needed beyond the cached HGNC
#' download. Symbols are normalized to their current approved form first,
#' since several sources report the symbol that was current when their
#' records were curated.
#'
#' \code{hgnc_id} is the identifier to merge on. \code{target_uniprot} and
#' \code{gene_symbol} are convenient but less reliable as keys: a handful
#' of genes carry several UniProt accessions and a handful of accessions
#' belong to more than one gene, and those cases are left \code{NA} rather
#' than resolved to an arbitrary choice.
#'
#' \code{compound_chembl_id} is filled for the sources that already
#' identify drugs by ChEMBL ID: ChEMBL itself, Open Targets, and DGIdb
#' records whose concept identifier is a ChEMBL one. TTD, the Broad Hub
#' and GtoPdb identify compounds by PubChem CID, InChIKey and IUPHAR
#' ligand ID, and are left \code{NA} here.
#'
#' @param results a named list of per-source \code{data.frame}s, as
#'   returned by \code{\link{queryDrugTargets}} or
#'   \code{\link{buildGenomeWideDrugTargetTable}}.
#' @param hgncTable an HGNC gene table as returned by
#'   \code{\link{getHgncGeneTable}}; downloaded (and cached) on demand when
#'   not supplied.
#' @param symbolMap an old-to-current symbol map as returned by
#'   \code{\link{buildHgncSymbolMap}}; built from \code{hgncTable} when not
#'   supplied.
#' @param verbose logical(1); if TRUE, report how many rows each source
#'   resolved.
#' @return \code{results} with the four columns above added to each
#'   table, in the same order and with the same number of rows.
#' @examples
#' \donttest{
#'   res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
#'                           sources = c("chembl", "opentargets", "dgidb"))
#'   res <- addCommonIds(res)
#'   merge(res$chembl, res$opentargets, by = c("hgnc_id", "compound_chembl_id"))
#' }
#' @seealso \code{\link{combineDrugTargets}}, \code{\link{getHgncGeneTable}}
#' @export
addCommonIds <- function(results, hgncTable = NULL, symbolMap = NULL,
                         verbose = FALSE) {
    if (!is.list(results) || is.data.frame(results))
        stop("'results' must be a named list of per-source data.frames, as ",
             "returned by queryDrugTargets() or ",
             "buildGenomeWideDrugTargetTable().")
    if (length(results) == 0L) return(results)
    unsupported <- setdiff(names(results), names(.dtiCommonIdSpec))
    if (length(unsupported))
        stop("addCommonIds() does not recognise source(s): ",
             paste(unsupported, collapse = ", "), ". Expected any of: ",
             paste(names(.dtiCommonIdSpec), collapse = ", "), ".")

    if (is.null(hgncTable)) hgncTable <- getHgncGeneTable()
    if (is.null(symbolMap)) symbolMap <- buildHgncSymbolMap(hgncTable)
    lk <- .dtiHgncLookups(hgncTable)

    for (src in names(results)) {
        df <- results[[src]]
        if (!is.data.frame(df) || nrow(df) == 0L) next
        spec <- .dtiCommonIdSpec[[src]]
        n <- nrow(df)
        col <- function(nm) if (!is.na(nm) && nm %in% names(df)) as.character(df[[nm]]) else NULL

        acc <- col(spec$uniprot)
        ens <- col(spec$ensembl)
        sym <- col(spec$symbol)
        if (!is.null(sym)) sym <- as.character(normalizeGeneSymbols(
            sym, symbolMap = symbolMap, hgncTable = hgncTable))

        ## Most specific identifier first: an accession or Ensembl id names
        ## one gene, a symbol is a label that may have been renamed.
        hgnc <- rep(NA_character_, n)
        fill <- function(hgnc, keys, lookup) {
            todo <- is.na(hgnc) & !is.na(keys)
            if (any(todo)) hgnc[todo] <- unname(lookup[keys[todo]])
            hgnc
        }
        if (!is.null(acc)) hgnc <- fill(hgnc, acc, lk$byAcc)
        if (!is.null(ens)) hgnc <- fill(hgnc, ens, lk$byEnsembl)
        if (!is.null(sym)) hgnc <- fill(hgnc, sym, lk$bySymbol)

        df$hgnc_id        <- hgnc
        df$gene_symbol    <- unname(lk$symbolOf[hgnc])
        df$target_uniprot <- if (!is.null(acc)) acc else unname(lk$accOf[hgnc])
        df$compound_chembl_id <- if (is.na(spec$compound)) rep(NA_character_, n)
                                 else .dtiAsChemblId(col(spec$compound))

        stopifnot(nrow(df) == n)  # lookups only - the grain must not change
        if (verbose)
            message("addCommonIds: ", src, " - ", sum(!is.na(df$hgnc_id)), "/", n,
                    " rows got an hgnc_id, ", sum(!is.na(df$compound_chembl_id)),
                    " a compound_chembl_id")
        results[[src]] <- df
    }
    results
}
