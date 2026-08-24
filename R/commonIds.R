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
#' compound entry names the source's own drug identifier column; whether
#' it actually holds a ChEMBL ID is decided per value by
#' \code{.dtiAsChemblId}, since DGIdb's is a ChEMBL id only when its
#' prefix says so.
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
#' A table that already carries an \code{hgnc_id} keeps it; only rows
#' where it is missing are resolved here. A genome-wide build records the
#' gene each row was queried from, and that is a better answer than
#' anything recoverable from the identifiers a source reports back.
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

        ## An hgnc_id that is already there is kept. buildGenomeWideDrugTargetTable()
        ## tags every row with the gene it queried *from*, which is more
        ## reliable than anything recoverable from what the source echoed
        ## back: re-deriving it recovers only a quarter of ChEMBL's rows,
        ## since many UniProt accessions name more than one gene, and for
        ## TTD it moves rows to the wrong gene outright (rows queried as
        ## ADRA1A resolve to HGNC:280 = ADRA1D, because ADRA1A is also a
        ## previous symbol of that gene). Gaps are still filled below.
        hgnc <- if ("hgnc_id" %in% names(df)) as.character(df$hgnc_id)
                else rep(NA_character_, n)
        kept <- sum(!is.na(hgnc))
        fill <- function(hgnc, keys, lookup) {
            todo <- is.na(hgnc) & !is.na(keys)
            if (any(todo)) hgnc[todo] <- unname(lookup[keys[todo]])
            hgnc
        }
        ## Most specific identifier first: an accession or Ensembl id names
        ## one gene, a symbol is a label that may have been renamed.
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
                    " rows got an hgnc_id", if (kept > 0L)
                        paste0(" (", kept, " kept as supplied)") else "",
                    ", ", sum(!is.na(df$compound_chembl_id)),
                    " a compound_chembl_id")
        results[[src]] <- df
    }
    results
}


## ---------------------------------------------------------------------
##  Horizontal (column-append) join
##
##  combineDrugTargets() stacks the sources; this puts them side by side,
##  one row per key with each source's columns appended. Every source
##  reports several rows per gene-drug pair - one per mechanism in ChEMBL,
##  per disease in Open Targets, per assay in GtoPdb - so a plain merge()
##  of six such tables multiplies those rows against each other. The many
##  values are therefore collapsed into one cell per key before anything is
##  joined, and the result is assembled by aligning each source against the
##  union of keys with match(), never by merge(), so the declared grain
##  holds by construction rather than by hope. .assertUniqueKey() then
##  checks it anyway.
## ---------------------------------------------------------------------

#' Collapse one source's rows to one row per key
#'
#' Returns a list with the unique key rows and, for each carried column, a
#' list of that key's distinct non-missing values.
#' @keywords internal
#' @noRd
.dtiCollapseBySource <- function(df, by, valueCols) {
    keyStr <- do.call(paste, c(unname(as.list(df[by])), list(sep = "\r")))
    idx <- split(seq_len(nrow(df)), factor(keyStr, levels = unique(keyStr)))
    keys <- df[vapply(idx, `[`, integer(1), 1L), by, drop = FALSE]
    vals <- lapply(valueCols, function(col) {
        x <- as.character(df[[col]])
        lapply(idx, function(i) unique(x[i][!is.na(x[i]) & nzchar(x[i])]))
    })
    names(vals) <- valueCols
    list(key = keyStr[vapply(idx, `[`, integer(1), 1L)], keys = keys, vals = vals)
}

#' Resolve a canonical column name to one source's own column name
#'
#' \code{columns} may name a shared column from the mapping table
#' (\code{"drug_name"}) or a source's own column (\code{"Drug_Name"}).
#' The first is expanded here into whatever that source calls it; the
#' second is passed through untouched, so calls written before the
#' mapping existed behave exactly as they did.
#' @keywords internal
#' @noRd
.dtiExpandColumns <- function(columns, colMap, src) {
    if (is.null(columns)) return(NULL)
    mapped <- colMap$column[colMap$canonical %in% columns & colMap$source == src]
    unique(c(columns, mapped[!is.na(mapped)]))
}

#' Join the per-source tables side by side, one row per gene-drug pair
#'
#' The horizontal counterpart to \code{\link{combineDrugTargets}}, which
#' stacks the sources on top of one another. This one places them side by
#' side: one row per key, with each source's own columns appended and
#' prefixed by the source name.
#'
#' Every source reports several rows for the same gene-drug pair - one per
#' mechanism in ChEMBL, one per disease in Open Targets, one per assay in
#' GtoPdb - so the values are collapsed to the distinct ones per key before
#' the sources are put together. Collapsed values are returned as
#' list-columns, which keep them addressable with \code{lengths()} and
#' \code{unlist()}; use \code{collapse = "string"} to get them as delimited
#' text instead, which is what you want for writing the table to a file.
#'
#' The key is set by \code{by}. The default, \code{c("hgnc_id",
#' "compound_chembl_id")}, gives one row per gene-drug pair; \code{by =
#' "hgnc_id"} gives one row per gene, with each source's drugs collapsed
#' into that gene's cell. Only three sources currently identify compounds
#' by ChEMBL ID (see \code{\link{addCommonIds}}), so a query including TTD,
#' the Broad Hub or GtoPdb keyed on a compound will drop their rows,
#' reporting how many; key on \code{"hgnc_id"} to keep them.
#'
#' \code{n_sources} counts the sources that returned a row for that key,
#' which is not the same as the sources that found a drug. A genome-wide
#' build queries every gene everywhere, so most of its genes carry a row
#' from every source with the drug columns empty. To count the sources
#' that actually found something, test the cells:
#' \code{lengths(x$chembl_Drug_Name) > 0}.
#'
#' @param results a named list of per-source \code{data.frame}s, as
#'   returned by \code{\link{queryDrugTargets}} or
#'   \code{\link{buildGenomeWideDrugTargetTable}}. Shared keys are added
#'   with \code{\link{addCommonIds}} if not already present.
#' @param by character vector of key columns, any of \code{"hgnc_id"},
#'   \code{"compound_chembl_id"}, \code{"gene_symbol"},
#'   \code{"target_uniprot"}.
#' @param columns character vector of columns to carry, given without the
#'   source prefix; the default carries every column each source has
#'   apart from the keys and \code{QueryIDs}. A shared column name from
#'   \code{colMap} (\code{"drug_name"}) selects whatever each source calls
#'   it, so you do not have to name all six variants; a source's own
#'   column name (\code{"Drug_Name"}) also works.
#' @param colMap the column mapping used to resolve shared names in
#'   \code{columns}; default \code{\link{drugTargetColumnMap}()}.
#' @param collapse \code{"list"} (default) for list-columns, or
#'   \code{"string"} to paste each cell's values together with \code{sep}.
#' @param sep separator used when \code{collapse = "string"}.
#' @param verbose logical(1); if TRUE, report per-source row counts and
#'   anything dropped for want of a key.
#' @return A \code{\link[S4Vectors]{DataFrame}} with one row per distinct
#'   \code{by} combination: the key columns first, then any gene identity
#'   the tables carry (\code{symbol}, \code{ensembl_gene_id} - once, not
#'   once per source), then \code{n_sources} and \code{sources}, then each
#'   source's carried columns prefixed with its name.
#' @examples
#' \donttest{
#'   res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
#'                           sources = c("chembl", "opentargets"))
#'   ## Each source's own column names...
#'   mergeDrugTargets(res, columns = c("Drug_Name", "Action_Type", "action_type"))
#'   ## ...or the shared names, which resolve to the same columns.
#'   mergeDrugTargets(res, columns = c("drug_name", "action"))
#' }
#' @seealso \code{\link{combineDrugTargets}} for the row-append form,
#'   \code{\link{addCommonIds}} for the keys themselves.
#' @export
mergeDrugTargets <- function(results, by = c("hgnc_id", "compound_chembl_id"),
                             columns = NULL, colMap = drugTargetColumnMap(),
                             collapse = c("list", "string"),
                             sep = " | ", verbose = FALSE) {
    collapse <- match.arg(collapse)
    by <- match.arg(by, .dtiCommonIdCols, several.ok = TRUE)
    if (!is.list(results) || is.data.frame(results))
        stop("'results' must be a named list of per-source data.frames, as ",
             "returned by queryDrugTargets() or ",
             "buildGenomeWideDrugTargetTable().")
    unsupported <- setdiff(names(results), names(.dtiCommonIdSpec))
    if (length(unsupported))
        stop("mergeDrugTargets() does not recognise source(s): ",
             paste(unsupported, collapse = ", "), ". Expected any of: ",
             paste(names(.dtiCommonIdSpec), collapse = ", "), ".")

    haveKeys <- vapply(results, function(d)
        is.data.frame(d) && all(.dtiCommonIdCols %in% names(d)), logical(1))
    if (length(haveKeys) && !all(haveKeys)) results <- addCommonIds(results)
    colMap <- .dtiValidateColumnMap(colMap)

    ## A genome-wide build tags every source's rows with the same gene
    ## identity, so prefixing it per source would repeat one answer once
    ## per source. It is collapsed like any other column but emitted once,
    ## beside the key.
    idCols <- setdiff(unique(unlist(lapply(results, function(d)
        if (is.data.frame(d)) intersect(.dtiIdentityCols, names(d)))),
        use.names = FALSE), c(.dtiCommonIdCols, "QueryIDs", by))

    parts <- list()
    for (src in names(results)) {
        df <- results[[src]]
        if (!is.data.frame(df) || nrow(df) == 0L) next
        keep <- Reduce(`&`, lapply(df[by], function(x) !is.na(x) & nzchar(x)))
        if (verbose && any(!keep))
            message("mergeDrugTargets: ", src, " - dropped ", sum(!keep), "/",
                    nrow(df), " row(s) with no ", paste(by, collapse = "/"))
        df <- df[keep, , drop = FALSE]
        if (nrow(df) == 0L) next
        valueCols <- setdiff(names(df), c(.dtiCommonIdCols, "QueryIDs", idCols))
        if (!is.null(columns))
            valueCols <- intersect(valueCols,
                                   .dtiExpandColumns(columns, colMap, src))
        parts[[src]] <- .dtiCollapseBySource(
            df, by, union(intersect(idCols, names(df)), valueCols))
        if (verbose)
            message("mergeDrugTargets: ", src, " - ", nrow(df), " row(s) -> ",
                    length(parts[[src]]$key), " key(s)")
    }
    if (length(parts) == 0L)
        return(S4Vectors::DataFrame(stats::setNames(
            replicate(length(by), character(0), simplify = FALSE), by)))

    ## Union of keys, in a deterministic order, then every source aligned
    ## onto it positionally - no merge(), so no fan-out is possible.
    keyRows <- unique(do.call(rbind, lapply(parts, `[[`, "keys")))
    ord <- do.call(order, unname(as.list(keyRows)))
    keyRows <- keyRows[ord, , drop = FALSE]
    keyAll <- do.call(paste, c(unname(as.list(keyRows)), list(sep = "\r")))
    n <- length(keyAll)

    out <- S4Vectors::DataFrame(keyRows, row.names = NULL)

    ## Gene identity: one column, holding what the sources agree on. Keyed
    ## on hgnc_id that is a single value per cell; keyed on a compound it
    ## is legitimately every gene that drug hits.
    for (col in idCols) {
        vs <- replicate(n, character(0), simplify = FALSE)
        for (p in parts) {
            if (is.null(p$vals[[col]])) next
            i <- match(keyAll, p$key)
            v <- p$vals[[col]][i]
            for (j in which(!is.na(i))) vs[[j]] <- union(vs[[j]], v[[j]])
        }
        out[[col]] <- I(vs)
    }

    present <- lapply(parts, function(p) !is.na(match(keyAll, p$key)))
    out$n_sources <- Reduce(`+`, lapply(present, as.integer))
    srcOf <- lapply(seq_len(n), function(i)
        names(parts)[vapply(present, `[`, logical(1), i)])
    out$sources <- I(srcOf)

    for (src in names(parts)) {
        p <- parts[[src]]
        i <- match(keyAll, p$key)
        for (col in setdiff(names(p$vals), idCols)) {
            v <- p$vals[[col]][i]
            v[is.na(i)] <- list(character(0))
            out[[paste0(src, "_", col)]] <- I(unname(v))
        }
    }

    if (identical(collapse, "string")) {
        for (col in names(out)) {
            if (!is.list(out[[col]])) next
            s <- vapply(out[[col]], function(x)
                if (length(x) == 0L) NA_character_ else paste(x, collapse = sep),
                character(1))
            out[[col]] <- s
        }
    }
    .assertUniqueKey(as.data.frame(keyRows), by, "mergeDrugTargets() result")
    out
}
