## =====================================================================
##  Cross-source mechanism-of-action (MOA) assembly
##
##  The sources this package supports disagree about MOA in two ways that
##  a naive column-append would silently paper over, so both are recorded
##  explicitly in the output rather than resolved away:
##
##  1. Grain. Open Targets attaches one mechanism statement that fans out
##     across several targets; ChEMBL/TTD/GtoPdb/DGIdb record it per
##     drug-target edge; the Broad Repurposing Hub's `moa` is a *drug*
##     level attribute - the identical " | "-joined string is repeated on
##     every one of that drug's target rows, so a particular MOA can never
##     be attributed to a particular target from Broad. Flattening those
##     together without a grain flag would invent target attributions that
##     the underlying data does not support, which is what `moa_scope`
##     exists to prevent.
##
##  2. Kind. ChEMBL/Open Targets/Broad give free-text mechanism prose
##     ("Bcr/Abl fusion protein inhibitor"); TTD/GtoPdb/DGIdb give a short
##     action vocabulary ("Inhibitor"). Both are legitimately "MOA" but
##     they are not interchangeable, hence `moa_kind`.
##
##  PubChem is absent by construction: it carries bioassay measurements
##  (an Active/Inactive call plus a potency value) and has no
##  mechanism-of-action concept at all - see listMoaSources().
##
##  Like combineDrugTargets(), assembleMoaTable() is a transformer over
##  queryDrugTargets()'s output rather than a querying function of its
##  own, which keeps all ID resolution and dispatch in one place and
##  makes the assembly logic a pure function that is testable offline.
## =====================================================================

#' Canonical column order of the long MOA table
#' @keywords internal
#' @noRd
.dtiMoaCols <- c("query_id", "drug_key", "drug_name", "source", "moa_text",
                 "moa_kind", "moa_scope", "action_type_raw", "action_type",
                 "target_symbol", "target_uniprot", "evidence_refs")

#' Per-source extraction spec for the long MOA table
#'
#' One entry per source, naming that source's own columns. \code{NA}
#' means the source has no equivalent column at all (not that it exists
#' but is empty). \code{split} names a column whose values pack several
#' MOA strings into one \code{" | "}-joined field (Broad only), which are
#' exploded to one row each. \code{scope} is the grain flag described in
#' this file's header; \code{"drug"} forces the target columns to
#' \code{NA} regardless of what the source's own target column holds.
#' @keywords internal
#' @noRd
.dtiMoaSourceSpec <- list(
    chembl = list(
        moa = "MOA", kind = "description", scope = "drug_target",
        action = "Action_Type", split = FALSE,
        symbol = NA, uniprot = "UniProt_ID", refs = NA,
        key = "chembl_id", name = "Drug_Name"),
    opentargets = list(
        moa = "mechanism_of_action", kind = "description", scope = "drug_target",
        action = "action_type", split = FALSE,
        symbol = "approved_symbol", uniprot = NA, refs = NA,
        key = "chembl_id", name = "drug_name"),
    broad = list(
        moa = "moa", kind = "description", scope = "drug",
        action = NA, split = TRUE,
        symbol = NA, uniprot = NA, refs = NA,
        key = "InChIKey", name = "pert_iname"),
    ttd = list(
        moa = "MOA", kind = "action_vocabulary", scope = "drug_target",
        action = "MOA", split = FALSE,
        symbol = "GeneName", uniprot = "Uniprot_acc", refs = NA,
        key = NA, name = "DrugName"),
    gtopdb = list(
        moa = "type", kind = "action_vocabulary", scope = "drug_target",
        action = "type", split = FALSE,
        symbol = "target_gene", uniprot = NA, refs = "refIds",
        key = NA, name = "ligandName"),
    dgidb = list(
        moa = "interaction_types", kind = "action_vocabulary", scope = "drug_target",
        action = "interaction_types", split = FALSE,
        symbol = "gene_name", uniprot = NA, refs = "sources",
        key = NA, name = "drug_name")
)

#' Action-term patterns, matched in order (first hit wins)
#'
#' The per-source action vocabularies are near-identical once case is
#' folded (ChEMBL \code{"INHIBITOR"}, TTD \code{"Inhibitor"}, DGIdb
#' \code{"inhibitor"}, GtoPdb \code{"Inhibitor"}/\code{"Inhibition"}), so
#' normalizing them is tractable in a way that normalizing free-text MOA
#' prose is not. Order is load-bearing rather than cosmetic:
#' \code{"antagonist"} contains \code{"agonist"} as a substring, and
#' \code{"inverse agonist"}/\code{"partial agonist"} are distinct terms
#' that must not collapse into plain \code{"agonist"} - so all three are
#' tested before it. Same priority-ordered keyword-bucketing approach as
#' \code{.ttdMoleculeType()}.
#' @keywords internal
#' @noRd
.dtiMoaActionPatterns <- c(
    "antagonist"      = "antagonist",
    "inverse agonist" = "inverse agonist",
    "partial agonist" = "partial agonist",
    "agonist"         = "agonist",
    "inhibit"         = "inhibitor",
    "block"           = "blocker",
    "activat"         = "activator",
    "modulat"         = "modulator",
    "stimulat"        = "stimulator",
    "substrate"       = "substrate",
    "chelat"          = "chelating agent",
    "disrupt"         = "disrupting agent",
    "releas"          = "releasing agent",
    "bind"            = "binder",
    "ligand"          = "ligand",
    "open"            = "opener"
)

#' Normalize a source's action term to a shared controlled vocabulary
#'
#' Returns \code{NA} for terms that match no pattern rather than guessing
#' - the raw term is always retained alongside in \code{action_type_raw},
#' so nothing is lost by declining to classify.
#' @keywords internal
#' @noRd
.dtiNormalizeActionType <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[x == ""] <- NA_character_
    out <- rep(NA_character_, length(x))
    for (pat in names(.dtiMoaActionPatterns)) {
        hit <- is.na(out) & !is.na(x) & grepl(pat, x, fixed = TRUE)
        out[hit] <- .dtiMoaActionPatterns[[pat]]
    }
    out
}

#' Cross-source drug key for one source's rows
#'
#' A ChEMBL molecule ID is the best available cross-source key (Open
#' Targets' drug IDs simply *are* ChEMBL IDs, and DGIdb's
#' \code{drug_aliases} carries one as a \code{"CHEMBL:CHEMBL348475"}
#' token), so it is preferred wherever obtainable; the Broad Hub has no
#' ChEMBL ID but does carry an InChIKey, which is a structure-level key
#' and the next best thing. TTD and GtoPdb carry neither, so they fall
#' back to a case-folded drug name. The fallback is deliberately weak and
#' documented as such: name matching across sources is unreliable, which
#' is exactly why \code{\link{assembleMoaTable}} does not merge rows on
#' \code{drug_key} - it is there for the caller to group on when they
#' judge it safe, not an identity claim made by this package.
#' @keywords internal
#' @noRd
.dtiMoaDrugKey <- function(df, src, spec) {
    keyCol <- spec$key
    key <- if (!is.na(keyCol) && keyCol %in% names(df)) {
        as.character(df[[keyCol]])
    } else if (src == "dgidb" && "drug_aliases" %in% names(df)) {
        ## "CHEMBL:CHEMBL348475; DRUGBANK:DB06240" -> "CHEMBL348475"
        vapply(seq_len(nrow(df)), function(i) {
            hit <- regmatches(df$drug_aliases[i],
                              regexpr("CHEMBL:CHEMBL[0-9]+", df$drug_aliases[i]))
            if (length(hit) == 1L) sub("^CHEMBL:", "", hit) else NA_character_
        }, character(1))
    } else {
        rep(NA_character_, nrow(df))
    }
    key[!is.na(key) & key == ""] <- NA_character_
    ## Name fallback, case-folded so "Imatinib"/"imatinib"/"IMATINIB" group.
    nm <- if (!is.na(spec$name) && spec$name %in% names(df))
        tolower(trimws(as.character(df[[spec$name]]))) else rep(NA_character_, nrow(df))
    ifelse(is.na(key), nm, key)
}

#' Which sources contribute mechanism-of-action data, and of what kind
#'
#' A quick-reference companion to \code{\link{assembleMoaTable}}: what
#' each source's MOA field actually contains, at what grain, and why
#' PubChem is excluded. All of it verified against live query output
#' rather than upstream documentation.
#'
#' @return A \code{data.frame} with one row per data source and columns
#'   \code{source}, \code{moa_column}, \code{moa_kind}, \code{moa_scope}
#'   and \code{note}.
#' @examples
#' listMoaSources()
#' @seealso \code{\link{assembleMoaTable}}
#' @export
listMoaSources <- function() {
    data.frame(
        source = c("chembl", "opentargets", "broad", "ttd", "gtopdb",
                   "dgidb", "pubchem"),
        moa_column = c("MOA", "mechanism_of_action", "moa", "MOA", "type",
                       "interaction_types", NA_character_),
        moa_kind = c("description", "description", "description",
                     "action_vocabulary", "action_vocabulary",
                     "action_vocabulary", NA_character_),
        moa_scope = c("drug_target", "drug_target", "drug", "drug_target",
                      "drug_target", "drug_target", NA_character_),
        note = c(
            "Free-text mechanism prose, per drug-target edge.",
            "Free-text prose; one statement can span several targets.",
            paste("Free-text prose but drug-level only - the same ' | '-joined",
                  "string repeats on every target row, so no per-target",
                  "attribution is possible."),
            paste("Short controlled vocabulary",
                  "(Inhibitor/Modulator/Antagonist/Agonist/...), ~5% missing."),
            paste("Controlled vocabulary; the paired 'action' column holds the",
                  "effect term (e.g. Inhibition)."),
            paste("Controlled vocabulary, frequently empty; DGIdb aggregates",
                  "other databases, so its rows may restate ChEMBL/TTD content."),
            paste("No mechanism-of-action data at all - bioassay measurements",
                  "only (Active/Inactive plus a potency value).")),
        stringsAsFactors = FALSE)
}

#' Assemble a long, grain-aware mechanism-of-action table
#'
#' Transforms \code{\link{queryDrugTargets}}'s per-source results into a
#' single long \code{data.frame} with one row per drug x source x MOA
#' statement x target, preserving each source's MOA string verbatim while
#' recording *what kind* of statement it is and *at what grain* it was
#' asserted. See this package's \code{listMoaSources()} for the per-source
#' summary.
#'
#' Two columns carry the distinctions that make the table safe to pool
#' across sources:
#' \describe{
#'   \item{\code{moa_scope}}{\code{"drug_target"} when the source asserted
#'     the mechanism for that specific drug-target pair, \code{"drug"}
#'     when it is a drug-level attribute that cannot be attributed to any
#'     one target (the Broad Repurposing Hub). \code{target_symbol} and
#'     \code{target_uniprot} are always \code{NA} on \code{"drug"} rows -
#'     filling them would fabricate an attribution the source does not
#'     make. \code{\link{moaByTarget}} drops those rows for you.}
#'   \item{\code{moa_kind}}{\code{"description"} for free-text mechanism
#'     prose, \code{"action_vocabulary"} for a short controlled term.
#'     Grouping or counting across the two without distinguishing them
#'     compares unlike things.}
#' }
#'
#' Rows are deduplicated but never merged across sources: the same
#' mechanism reported by ChEMBL and by Open Targets stays as two rows,
#' since collapsing them would require an identity claim this package
#' deliberately does not make (see \code{.dtiMoaDrugKey}). Free-text MOA
#' strings are likewise left verbatim - real synonymy resolution needs an
#' ontology, so only \code{action_type} is normalized (see
#' \code{.dtiNormalizeActionType}).
#'
#' Sources also disagree on how they name a target: ChEMBL's REST output
#' is UniProt-accession-keyed with no gene symbol at all, while most
#' others carry only a symbol, so \code{target_symbol} alone will not line
#' ChEMBL up with the rest. \code{resolveGeneSymbol = TRUE} fills the
#' missing symbols from \code{target_uniprot} in one batched lookup - the
#' same opt-in, one-extra-round-trip treatment
#' \code{\link{combineDrugTargets}} gives its \code{gene_symbol} column,
#' and off by default for the same reason (assembly stays network-free
#' unless asked otherwise).
#'
#' @param results a named list as returned by \code{\link{queryDrugTargets}}
#'   (or any similarly-shaped named list of per-source data.frames -
#'   \code{query_id} falls back to each source's own \code{QueryIDs}
#'   column when \code{results} has no \code{"resolved"} attribute).
#'   Sources with no MOA data, and sources present but matching zero rows,
#'   simply contribute nothing.
#' @param resolveGeneSymbol logical(1); if \code{TRUE}, fill
#'   otherwise-\code{NA} \code{target_symbol} values by resolving
#'   \code{target_uniprot} via \code{\link{getUniprotMapping}} (one extra
#'   network round trip; default \code{FALSE}). Rows whose accession does
#'   not resolve keep \code{NA}.
#' @param taxId integer(1) passed to \code{\link{getUniprotMapping}} when
#'   \code{resolveGeneSymbol = TRUE} (default 9606L = human).
#' @return A \code{data.frame} with columns \code{query_id},
#'   \code{drug_key}, \code{drug_name}, \code{source}, \code{moa_text},
#'   \code{moa_kind}, \code{moa_scope}, \code{action_type_raw},
#'   \code{action_type}, \code{target_symbol}, \code{target_uniprot} and
#'   \code{evidence_refs}.
#' @examples
#' ## A pure transformer, so it runs on any correctly-shaped per-source
#' ## frames - no network needed. Note Broad's drug-level MOA explodes to
#' ## one row per mechanism and keeps no target attribution.
#' chembl <- data.frame(QueryIDs = "CHEMBL941", chembl_id = "CHEMBL941",
#'                      Drug_Name = "IMATINIB",
#'                      MOA = "Tyrosine-protein kinase ABL inhibitor",
#'                      Action_Type = "INHIBITOR", UniProt_ID = "P00519")
#' broad <- data.frame(QueryIDs = "imatinib", pert_iname = "imatinib",
#'                     InChIKey = "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
#'                     moa = "Bcr-Abl kinase inhibitor | KIT inhibitor",
#'                     target_gene = "ABL1")
#' moa <- assembleMoaTable(list(chembl = chembl, broad = broad))
#' moa[, c("source", "moa_text", "moa_kind", "moa_scope", "target_symbol")]
#' \donttest{
#'   ## the usual path: straight from queryDrugTargets()
#'   res <- queryDrugTargets(list(molType = "cmp", idType = "chembl_id",
#'                                ids = "CHEMBL941"), sources = "chembl")
#'   assembleMoaTable(res, resolveGeneSymbol = TRUE)
#' }
#' @seealso \code{\link{queryMoa}}, \code{\link{moaWide}},
#'   \code{\link{moaByTarget}}, \code{\link{listMoaSources}},
#'   \code{\link{combineDrugTargets}}
#' @export
assembleMoaTable <- function(results, resolveGeneSymbol = FALSE, taxId = 9606L) {
    emptyOut <- as.data.frame(stats::setNames(
        replicate(length(.dtiMoaCols), character(0), simplify = FALSE), .dtiMoaCols))
    if (length(results) == 0L) return(emptyOut)

    unsupported <- setdiff(names(results), names(.dtiMoaSourceSpec))
    if (length(unsupported))
        stop("assembleMoaTable() has no MOA mapping for source(s): ",
             paste(unsupported, collapse = ", "), ". Sources carrying MOA data ",
             "are: ", paste(names(.dtiMoaSourceSpec), collapse = ", "),
             " (see listMoaSources() for why others are excluded).")

    resolvedAttr <- attr(results, "resolved")

    rows <- lapply(names(results), function(src) {
        df <- results[[src]]
        spec <- .dtiMoaSourceSpec[[src]]
        if (is.null(df) || nrow(df) == 0L) return(NULL)  ## resolved but matched
                                                          ## nothing; contributes
                                                          ## no rows rather than
                                                          ## erroring
        if (!spec$moa %in% names(df))
            stop("assembleMoaTable(): expected column '", spec$moa, "' not found ",
                 "in results$", src, " - is this really a ", src, " result from ",
                 "queryDrugTargets()?")

        pick <- function(col) if (!is.na(col) && col %in% names(df))
            as.character(df[[col]]) else rep(NA_character_, nrow(df))

        queryId <- if (!is.null(resolvedAttr[[src]])) {
            orig <- names(resolvedAttr[[src]])
            names(orig) <- unname(resolvedAttr[[src]])
            unname(orig[df$QueryIDs])
        } else if ("QueryIDs" %in% names(df)) {
            as.character(df$QueryIDs)
        } else rep(NA_character_, nrow(df))

        out <- data.frame(
            query_id        = queryId,
            drug_key        = .dtiMoaDrugKey(df, src, spec),
            drug_name       = pick(spec$name),
            source          = src,
            moa_text        = pick(spec$moa),
            moa_kind        = spec$kind,
            moa_scope       = spec$scope,
            action_type_raw = pick(spec$action),
            ## A drug-level MOA cannot be pinned to a target - see the
            ## file header. Blanking these here (rather than trusting each
            ## caller to notice moa_scope) is what keeps the invariant.
            target_symbol   = if (spec$scope == "drug") NA_character_ else pick(spec$symbol),
            target_uniprot  = if (spec$scope == "drug") NA_character_ else pick(spec$uniprot),
            evidence_refs   = pick(spec$refs),
            stringsAsFactors = FALSE)

        ## Broad packs several mechanisms into one " | "-joined field; one
        ## row each, mirroring .brhExplodeTargets()'s handling of that
        ## source's equally-packed target column.
        if (isTRUE(spec$split)) {
            parts <- strsplit(ifelse(is.na(out$moa_text), "", out$moa_text),
                              "|", fixed = TRUE)
            n <- lengths(parts)
            n[n == 0L] <- 1L      ## an empty MOA still yields its one NA row
            out <- out[rep(seq_len(nrow(out)), n), , drop = FALSE]
            flat <- trimws(unlist(lapply(parts, function(p)
                if (length(p) == 0L) NA_character_ else p), use.names = FALSE))
            out$moa_text <- flat
        }

        out$moa_text[!is.na(out$moa_text) & out$moa_text == ""] <- NA_character_
        out$action_type_raw[!is.na(out$action_type_raw) &
                            out$action_type_raw == ""] <- NA_character_
        out$action_type <- .dtiNormalizeActionType(out$action_type_raw)
        ## Drop rows carrying no mechanism information whatsoever - a
        ## target hit with a blank MOA belongs in queryDrugTargets()'s
        ## output, not in a MOA table.
        out <- out[!(is.na(out$moa_text) & is.na(out$action_type_raw)), , drop = FALSE]
        if (nrow(out) == 0L) return(NULL)
        out[, .dtiMoaCols, drop = FALSE]
    })

    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) return(emptyOut)
    out <- unique(do.call(rbind, rows))
    rownames(out) <- NULL

    if (isTRUE(resolveGeneSymbol)) {
        ## Source-agnostic on purpose: any row that names its target only
        ## by accession gets a symbol, which today is ChEMBL's whole output
        ## plus TTD's family-level rows, but needs no per-source special
        ## casing if another source later does the same.
        fill <- is.na(out$target_symbol) & !is.na(out$target_uniprot)
        accs <- unique(out$target_uniprot[fill])
        if (length(accs) > 0L) {
            sym <- .resolveGeneIds(accs, idType = "uniprot", to = "symbol",
                                   taxId = taxId)
            out$target_symbol[fill] <- unname(sym[out$target_uniprot[fill]])
        }
    }
    out
}

#' Query one or more sources and assemble their MOA table in one call
#'
#' Convenience wrapper: \code{\link{queryDrugTargets}} followed by
#' \code{\link{assembleMoaTable}}. Use the two separately when you also
#' want each source's full native columns, which this wrapper discards.
#'
#' @param queryBy a \code{list(molType, idType, ids)} as accepted by
#'   \code{\link{queryDrugTargets}}.
#' @param sources character vector of sources to query; defaults to every
#'   source that carries MOA data (see \code{\link{listMoaSources}}).
#' @param resolveGeneSymbol logical(1); passed to
#'   \code{\link{assembleMoaTable}}.
#' @param taxId integer(1); passed to \code{\link{assembleMoaTable}}.
#' @param ... further arguments passed to \code{\link{queryDrugTargets}},
#'   e.g. \code{ttdDbPath}, \code{brhDbPath}, \code{gtoPdbDbPath},
#'   \code{unichemDbPath}, \code{verbose}.
#' @return The \code{data.frame} described in \code{\link{assembleMoaTable}}.
#' @note Querying compounds by \code{idType = "name"} routes through
#'   \code{\link{queryDrugTargets}}'s compound resolution, which reaches
#'   ChEMBL and Open Targets only via a ChEMBL ID - so those two sources
#'   need \code{unichemDbPath} (see \code{\link{buildUnichemDb}}) when the
#'   query starts from a name. Without it they are skipped and the table
#'   is assembled from the remaining sources.
#' @examples
#' \donttest{
#'   queryMoa(list(molType = "cmp", idType = "chembl_id", ids = "CHEMBL941"),
#'            sources = c("chembl", "opentargets"))
#' }
#' @seealso \code{\link{assembleMoaTable}}
#' @export
queryMoa <- function(queryBy, sources = names(.dtiMoaSourceSpec),
                     resolveGeneSymbol = FALSE, taxId = 9606L, ...) {
    assembleMoaTable(queryDrugTargets(queryBy, sources = sources, ...),
                     resolveGeneSymbol = resolveGeneSymbol, taxId = taxId)
}

#' Collapse a long MOA table to one row per drug
#'
#' Reshapes \code{\link{assembleMoaTable}}'s output column-wise, with
#' list-columns where a drug has several distinct values. Both MOA kinds
#' are kept but in separate columns, since free-text prose and controlled
#' action terms are not interchangeable (see \code{\link{assembleMoaTable}}).
#'
#' @param x a \code{data.frame} as returned by \code{\link{assembleMoaTable}}.
#' @return A \code{data.frame} with one row per \code{query_id} x
#'   \code{drug_key}, columns \code{query_id}, \code{drug_key},
#'   \code{drug_name}, \code{n_sources}, and the list-columns
#'   \code{sources}, \code{moa_description}, \code{moa_action_vocabulary}
#'   and \code{action_type}.
#' @examples
#' chembl <- data.frame(QueryIDs = "CHEMBL941", chembl_id = "CHEMBL941",
#'                      Drug_Name = "IMATINIB",
#'                      MOA = "Tyrosine-protein kinase ABL inhibitor",
#'                      Action_Type = "INHIBITOR", UniProt_ID = "P00519")
#' wide <- moaWide(assembleMoaTable(list(chembl = chembl)))
#' wide[, c("drug_key", "drug_name", "n_sources")]
#' wide$moa_description[[1]]
#' @seealso \code{\link{assembleMoaTable}}, \code{\link{moaByTarget}}
#' @export
moaWide <- function(x) {
    .dtiAssertMoaTable(x)
    if (nrow(x) == 0L) {
        out <- data.frame(query_id = character(0), drug_key = character(0),
                          drug_name = character(0), n_sources = integer(0),
                          stringsAsFactors = FALSE)
        out$sources <- list()
        out$moa_description <- list()
        out$moa_action_vocabulary <- list()
        out$action_type <- list()
        return(out)
    }
    grp <- paste(x$query_id, x$drug_key, sep = "\r")
    idx <- split(seq_len(nrow(x)), factor(grp, levels = unique(grp)))
    uniqNoNA <- function(v) {
        v <- unique(v[!is.na(v)])
        if (length(v) == 0L) character(0) else v
    }
    out <- data.frame(
        query_id  = vapply(idx, function(i) x$query_id[i][1L], character(1)),
        drug_key  = vapply(idx, function(i) x$drug_key[i][1L], character(1)),
        drug_name = vapply(idx, function(i) {
            nm <- uniqNoNA(x$drug_name[i]); if (length(nm)) nm[1L] else NA_character_
        }, character(1)),
        n_sources = vapply(idx, function(i) length(unique(x$source[i])), integer(1)),
        stringsAsFactors = FALSE)
    out$sources <- unname(lapply(idx, function(i) unique(x$source[i])))
    out$moa_description <- unname(lapply(idx, function(i)
        uniqNoNA(x$moa_text[i][x$moa_kind[i] == "description"])))
    out$moa_action_vocabulary <- unname(lapply(idx, function(i)
        uniqNoNA(x$moa_text[i][x$moa_kind[i] == "action_vocabulary"])))
    out$action_type <- unname(lapply(idx, function(i) uniqNoNA(x$action_type[i])))
    rownames(out) <- NULL
    out
}

#' Keep only MOA rows that are attributable to a specific target
#'
#' Drops \code{moa_scope == "drug"} rows - currently the Broad
#' Repurposing Hub's, whose MOA is a drug-level attribute that cannot be
#' pinned to any one target (see \code{\link{assembleMoaTable}}). Use this
#' before any per-target grouping, so drug-level statements are not
#' silently counted as target-level evidence.
#'
#' @param x a \code{data.frame} as returned by \code{\link{assembleMoaTable}}.
#' @param verbose logical(1); if \code{TRUE} (default) report how many
#'   rows were dropped and from which sources.
#' @return \code{x} with only \code{moa_scope == "drug_target"} rows.
#' @examples
#' chembl <- data.frame(QueryIDs = "CHEMBL941", chembl_id = "CHEMBL941",
#'                      Drug_Name = "IMATINIB",
#'                      MOA = "Tyrosine-protein kinase ABL inhibitor",
#'                      Action_Type = "INHIBITOR", UniProt_ID = "P00519")
#' broad <- data.frame(QueryIDs = "imatinib", pert_iname = "imatinib",
#'                     InChIKey = "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
#'                     moa = "Bcr-Abl kinase inhibitor", target_gene = "ABL1")
#' moa <- assembleMoaTable(list(chembl = chembl, broad = broad))
#' ## the Broad row is dropped - its MOA is not attributable to a target
#' moaByTarget(moa)[, c("source", "moa_text", "target_uniprot")]
#' @seealso \code{\link{assembleMoaTable}}, \code{\link{moaWide}}
#' @export
moaByTarget <- function(x, verbose = TRUE) {
    .dtiAssertMoaTable(x)
    drop <- x$moa_scope == "drug"
    if (isTRUE(verbose) && any(drop))
        message("moaByTarget(): dropped ", sum(drop), " drug-level MOA row(s) from ",
                paste(sort(unique(x$source[drop])), collapse = ", "),
                " - those sources do not attribute a mechanism to a specific target.")
    out <- x[!drop, , drop = FALSE]
    rownames(out) <- NULL
    out
}

#' Validate that an object looks like assembleMoaTable() output
#' @keywords internal
#' @noRd
.dtiAssertMoaTable <- function(x) {
    if (!is.data.frame(x))
        stop("Expected a data.frame as returned by assembleMoaTable(), got ",
             class(x)[1L], ".")
    missing <- setdiff(.dtiMoaCols, names(x))
    if (length(missing))
        stop("Not an assembleMoaTable() result - missing column(s): ",
             paste(missing, collapse = ", "), ".")
    invisible(TRUE)
}
