## =====================================================================
##  Cross-source mechanism-of-action (MOA) assembly
##
##  WHAT A MOA IS HERE
##
##  A mechanism of action describes the molecular/biochemical action of a
##  *drug* - what it does and, almost always, to which molecular target
##  ("Bcr/Abl fusion protein inhibitor"). It is therefore a property of
##  the drug, not of a drug-target edge. Two measurements on ChEMBL's
##  drug_mechanism table pin this down:
##
##    * 48% of distinct MOA terms are used by more than one drug
##      ("Histamine H1 receptor antagonist" covers 55 of them), so a MOA
##      is a shared functional descriptor, close to an ontology term.
##    * a MOA term is nonetheless usually tied to one target: only 26 of
##      5,641 (drug, MOA) pairs span more than one target (0.5%).
##
##  So the natural shape is two relations, not one flat table: the drug's
##  MOA terms (assembleMoaTable), and the separate link from a MOA to the
##  target(s) a source says it acts through (assembleMoaTargets). A
##  source that states a MOA without naming a target - the Broad
##  Repurposing Hub - is then simply a source that contributes no link
##  rows, rather than a special case needing a flag.
##
##  WHICH SOURCES ACTUALLY CARRY MOA
##
##  Only ChEMBL, Open Targets and the Broad Hub. TTD, GtoPdb and DGIdb
##  label a column "MOA"/"type"/"interaction_types" but store an *action
##  type* - the verb alone, with no target named. The cardinalities are
##  decisive, and ChEMBL settles the question at schema level by carrying
##  both concepts as separate fields:
##
##    ChEMBL mechanism_of_action   1,589 distinct   <- MOA
##    ChEMBL action_type              33 distinct   <- action type
##    TTD    MOA                      48 distinct   <- action type
##    GtoPdb type                     10 distinct   <- action type (and
##                                                    partly modality:
##                                                    "Antibody", "None")
##
##  Those three sources are therefore excluded here rather than pooled
##  under a MOA heading, where their terms would swamp any grouping or
##  frequency count ("Inhibitor" would be the commonest "MOA" in the
##  table). Nothing is lost: combineDrugTargets()'s `action` column
##  already carries that dimension for all six sources, which is where it
##  belongs. PubChem has neither - bioassay measurements only.
##
##  Mode of Action (MoA - the broader cellular/organism-level consequence,
##  as distinct from MOA) is not modelled: no source carries it
##  systematically. 94% of ChEMBL's MOA strings are [target] + [action],
##  and the rest are still molecular or modality ("Iron chelating agent",
##  "Surfactant agent"). The only MoA-flavoured content anywhere is
##  ChEMBL's free-text mechanism_comment ("Role in regulating gastric
##  secretion, M3 likely involved"), which is carried through as an
##  optional passthrough column rather than modelled as a concept.
##
##  Like combineDrugTargets(), these are transformers over
##  queryDrugTargets()'s output rather than querying functions of their
##  own, which keeps ID resolution in one place and makes the assembly
##  logic pure and testable offline.
## =====================================================================

#' Canonical column order of the drug -> MOA table
#' @keywords internal
#' @noRd
.dtiMoaCols <- c("query_id", "drug_key", "drug_name", "source", "moa_text",
                 "action_type_raw", "action_type")

#' Canonical column order of the MOA -> target link table
#' @keywords internal
#' @noRd
.dtiMoaTargetCols <- c("query_id", "drug_key", "source", "moa_text",
                       "target_symbol", "target_uniprot", "mechanism_comment")

#' Per-source extraction spec
#'
#' Only the three sources that carry genuine MOA terms appear here; see
#' this file's header for why TTD/GtoPdb/DGIdb do not. \code{NA} means
#' the source has no equivalent column at all. \code{split} names a
#' source whose MOA field packs several terms into one \code{" | "}
#' -joined string (Broad only). \code{linksTargets} says whether the
#' source ties its MOA to a named target - \code{FALSE} for Broad, whose
#' \code{moa} is stated for the drug with no per-target attribution, so
#' it contributes to the MOA table but not the link table.
#' @keywords internal
#' @noRd
.dtiMoaSourceSpec <- list(
    chembl = list(
        moa = "MOA", action = "Action_Type", split = FALSE, linksTargets = TRUE,
        symbol = NA, uniprot = "UniProt_ID",
        comment = "mechanism.mechanism_comment",
        key = "chembl_id", name = "Drug_Name"),
    opentargets = list(
        moa = "mechanism_of_action", action = "action_type", split = FALSE,
        linksTargets = TRUE,
        symbol = "approved_symbol", uniprot = NA, comment = NA,
        key = "chembl_id", name = "drug_name"),
    broad = list(
        moa = "moa", action = NA, split = TRUE, linksTargets = FALSE,
        symbol = NA, uniprot = NA, comment = NA,
        key = "InChIKey", name = "pert_iname")
)

#' Action-term patterns, matched in order (first hit wins)
#'
#' Normalizes ChEMBL's and Open Targets' action vocabularies to shared
#' terms. Order is load-bearing rather than cosmetic:
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
    "cross-link"      = "cross-linking agent",
    "bind"            = "binder",
    "ligand"          = "ligand",
    "open"            = "opener"
)

#' Normalize an action term to a shared controlled vocabulary
#'
#' Returns \code{NA} for terms matching no pattern rather than guessing -
#' the source's own term is always kept in \code{action_type_raw}, so
#' declining to classify loses nothing.
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
#' Targets' drug IDs simply *are* ChEMBL IDs), so it is preferred; the
#' Broad Hub has no ChEMBL ID but does carry an InChIKey, a
#' structure-level key and the next best thing. The case-folded-name
#' fallback is deliberately weak, and documented as such: name matching
#' across sources is unreliable, which is why neither assembler merges
#' rows on \code{drug_key}. It is there for the caller to group on when
#' they judge it safe, not an identity claim made by this package.
#' @keywords internal
#' @noRd
.dtiMoaDrugKey <- function(df, spec) {
    key <- if (!is.na(spec$key) && spec$key %in% names(df)) {
        as.character(df[[spec$key]])
    } else rep(NA_character_, nrow(df))
    key[!is.na(key) & key == ""] <- NA_character_
    nm <- if (!is.na(spec$name) && spec$name %in% names(df))
        tolower(trimws(as.character(df[[spec$name]]))) else rep(NA_character_, nrow(df))
    ifelse(is.na(key), nm, key)
}

#' Shared per-source extraction used by both assemblers
#'
#' Returns one row per source row x MOA term, with every column both
#' assemblers need; each public function then selects and dedupes its
#' own. Keeping this in one place is what guarantees the two tables
#' cannot disagree about which MOA terms exist.
#' @keywords internal
#' @noRd
.dtiMoaExtract <- function(results) {
    unsupported <- setdiff(names(results), names(.dtiMoaSourceSpec))
    if (length(unsupported))
        stop("MOA assembly does not use source(s): ",
             paste(unsupported, collapse = ", "),
             ". Only ChEMBL, Open Targets and the Broad Repurposing Hub carry ",
             "mechanism-of-action terms; see listMoaSources() for what the ",
             "others store instead, and combineDrugTargets() for their ",
             "action-type data.")

    resolvedAttr <- attr(results, "resolved")

    rows <- lapply(names(results), function(src) {
        df <- results[[src]]
        spec <- .dtiMoaSourceSpec[[src]]
        if (is.null(df) || nrow(df) == 0L) return(NULL)  ## resolved but matched
                                                          ## nothing: contributes
                                                          ## no rows, no error
        if (!spec$moa %in% names(df))
            stop("MOA assembly: expected column '", spec$moa, "' not found in ",
                 "results$", src, " - is this really a ", src, " result from ",
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
            query_id          = queryId,
            drug_key          = .dtiMoaDrugKey(df, spec),
            drug_name         = pick(spec$name),
            source            = src,
            moa_text          = pick(spec$moa),
            action_type_raw   = pick(spec$action),
            target_symbol     = pick(spec$symbol),
            target_uniprot    = pick(spec$uniprot),
            mechanism_comment = pick(spec$comment),
            stringsAsFactors  = FALSE)

        ## Broad packs several MOA terms into one " | "-joined field; one
        ## row each, mirroring .brhExplodeTargets()'s handling of that
        ## source's equally-packed target column.
        if (isTRUE(spec$split)) {
            parts <- strsplit(ifelse(is.na(out$moa_text), "", out$moa_text),
                              "|", fixed = TRUE)
            n <- lengths(parts)
            n[n == 0L] <- 1L      ## an empty MOA still yields its one NA row
            out <- out[rep(seq_len(nrow(out)), n), , drop = FALSE]
            out$moa_text <- trimws(unlist(lapply(parts, function(p)
                if (length(p) == 0L) NA_character_ else p), use.names = FALSE))
        }

        ## A source that does not tie its MOA to a target must never
        ## appear to: blanking here is what keeps the link table honest.
        if (!isTRUE(spec$linksTargets)) {
            out$target_symbol <- NA_character_
            out$target_uniprot <- NA_character_
        }

        blank <- function(v) { v[!is.na(v) & v == ""] <- NA_character_; v }
        out$moa_text <- blank(out$moa_text)
        out$action_type_raw <- blank(out$action_type_raw)
        out$mechanism_comment <- blank(out$mechanism_comment)
        ## A row with no MOA term is not MOA data - a bare target hit
        ## belongs in queryDrugTargets()'s output, not here.
        out <- out[!is.na(out$moa_text), , drop = FALSE]
        if (nrow(out) == 0L) return(NULL)
        out$action_type <- .dtiNormalizeActionType(out$action_type_raw)
        out
    })

    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) return(NULL)
    do.call(rbind, rows)
}

#' Which sources carry mechanism-of-action terms, and what the others store
#'
#' Reference table behind \code{\link{assembleMoaTable}}'s source
#' selection, including the distinct-value counts that separate a MOA
#' field from an action-type field. All figures measured from the sources
#' themselves rather than taken from their documentation.
#'
#' @return A \code{data.frame} with one row per data source and columns
#'   \code{source}, \code{column}, \code{kind}, \code{n_distinct},
#'   \code{used_for_moa} and \code{note}.
#' @examples
#' listMoaSources()
#' @seealso \code{\link{assembleMoaTable}}, \code{\link{combineDrugTargets}}
#' @export
listMoaSources <- function() {
    data.frame(
        source = c("chembl", "opentargets", "broad", "ttd", "gtopdb",
                   "dgidb", "pubchem"),
        column = c("MOA", "mechanism_of_action", "moa", "MOA", "type",
                   "interaction_types", NA_character_),
        kind = c("mechanism of action", "mechanism of action",
                 "mechanism of action", "action type", "action type",
                 "action type", "none"),
        n_distinct = c(1589L, NA_integer_, 1612L, 48L, 10L, NA_integer_,
                       NA_integer_),
        used_for_moa = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
        note = c(
            paste("Reference model: names a molecular target plus an action",
                  "(94% of terms). Carries action_type as a separate field,",
                  "which is what makes the two concepts distinguishable."),
            paste("Same shape as ChEMBL and largely ChEMBL-derived; counted",
                  "per query rather than in bulk, hence no n_distinct."),
            paste("Genuine MOA terms, but stated for the drug with no",
                  "per-target attribution, so it contributes no rows to",
                  "assembleMoaTargets()."),
            paste("Labelled MOA but stores an action type - 48 values such as",
                  "Inhibitor/Agonist/Blocker, none naming a target. Use",
                  "combineDrugTargets() for this dimension."),
            paste("Action type, and partly drug modality ('Antibody',",
                  "'Fusion protein', 'None'). Use combineDrugTargets()."),
            paste("Interaction-type vocabulary, frequently empty, and DGIdb",
                  "aggregates other databases. Use combineDrugTargets()."),
            paste("No mechanism data at all - bioassay measurements only",
                  "(Active/Inactive plus a potency value).")),
        stringsAsFactors = FALSE)
}

#' Assemble a drug -> mechanism-of-action table
#'
#' Transforms \code{\link{queryDrugTargets}}'s per-source results into
#' the MOA terms each source states for each drug, one row per drug x
#' source x MOA term. Because a mechanism of action is a property of the
#' *drug*, the target(s) a source says the mechanism acts through are a
#' separate relation - see \code{\link{assembleMoaTargets}}.
#'
#' Only ChEMBL, Open Targets and the Broad Repurposing Hub are used:
#' TTD, GtoPdb and DGIdb label a column "MOA"/"type"/"interaction_types"
#' but store an *action type* (the verb alone, no target named), and
#' pooling those under a MOA heading would let "Inhibitor" dominate any
#' grouping or frequency count. \code{\link{listMoaSources}} gives the
#' evidence; \code{\link{combineDrugTargets}}'s \code{action} column
#' already carries that dimension for all six sources.
#'
#' MOA strings are kept verbatim - genuine synonymy across free text
#' would need an ontology. Only \code{action_type} is normalized, and
#' only for the sources that state one, since those vocabularies are
#' near-identical once case is folded.
#'
#' Rows are deduplicated but never merged across sources: the same
#' mechanism reported by ChEMBL and by Open Targets stays as two rows,
#' because collapsing them would require a cross-source identity claim
#' this package deliberately does not make.
#'
#' @param results a named list as returned by \code{\link{queryDrugTargets}}
#'   (or any similarly-shaped named list of per-source data.frames -
#'   \code{query_id} falls back to each source's own \code{QueryIDs}
#'   column when \code{results} has no \code{"resolved"} attribute).
#'   Sources present but matching zero rows contribute nothing; sources
#'   that carry no MOA are rejected with a pointer to
#'   \code{\link{listMoaSources}}.
#' @return A \code{data.frame} with columns \code{query_id},
#'   \code{drug_key}, \code{drug_name}, \code{source}, \code{moa_text},
#'   \code{action_type_raw} and \code{action_type}.
#' @examples
#' ## A pure transformer, so it runs on any correctly-shaped per-source
#' ## frames - no network needed. ChEMBL's two records for one mechanism
#' ## against two targets collapse to a single MOA row.
#' chembl <- data.frame(
#'     QueryIDs = c("CHEMBL941", "CHEMBL941"), chembl_id = "CHEMBL941",
#'     Drug_Name = "IMATINIB", MOA = "Bcr/Abl fusion protein inhibitor",
#'     Action_Type = "INHIBITOR", UniProt_ID = c("P00519", "P11274"))
#' broad <- data.frame(
#'     QueryIDs = "imatinib", pert_iname = "imatinib",
#'     InChIKey = "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
#'     moa = "Bcr-Abl kinase inhibitor | KIT inhibitor", target_gene = "ABL1")
#' assembleMoaTable(list(chembl = chembl, broad = broad))
#' @seealso \code{\link{assembleMoaTargets}}, \code{\link{queryMoa}},
#'   \code{\link{moaWide}}, \code{\link{listMoaSources}},
#'   \code{\link{combineDrugTargets}}
#' @export
assembleMoaTable <- function(results) {
    emptyOut <- as.data.frame(stats::setNames(
        replicate(length(.dtiMoaCols), character(0), simplify = FALSE),
        .dtiMoaCols))
    if (length(results) == 0L) return(emptyOut)
    ext <- .dtiMoaExtract(results)
    if (is.null(ext)) return(emptyOut)
    out <- unique(ext[, .dtiMoaCols, drop = FALSE])
    rownames(out) <- NULL
    out
}

#' Assemble the link from a MOA term to the target(s) it acts through
#'
#' The companion relation to \code{\link{assembleMoaTable}}: one row per
#' drug x source x MOA term x target, for those sources that tie a
#' mechanism to a named target. The Broad Repurposing Hub states its MOA
#' for the drug with no per-target attribution and so contributes no rows
#' here - which is the point of keeping the two relations separate rather
#' than flagging drug-level rows inside one table.
#'
#' Sources also disagree on how they name a target: ChEMBL's REST output
#' is UniProt-accession-keyed with no gene symbol at all, while Open
#' Targets carries only a symbol. \code{resolveGeneSymbol = TRUE} fills
#' the missing symbols from \code{target_uniprot} in one batched lookup -
#' the same opt-in, one-extra-round-trip treatment
#' \code{\link{combineDrugTargets}} gives its \code{gene_symbol} column,
#' and off by default for the same reason.
#'
#' @param results a named list as returned by \code{\link{queryDrugTargets}};
#'   see \code{\link{assembleMoaTable}}.
#' @param resolveGeneSymbol logical(1); if \code{TRUE}, fill
#'   otherwise-\code{NA} \code{target_symbol} values by resolving
#'   \code{target_uniprot} via \code{\link{getUniprotMapping}} (one extra
#'   network round trip; default \code{FALSE}). Accessions that do not
#'   resolve keep \code{NA}.
#' @param taxId integer(1) passed to \code{\link{getUniprotMapping}} when
#'   \code{resolveGeneSymbol = TRUE} (default 9606L = human).
#' @return A \code{data.frame} with columns \code{query_id},
#'   \code{drug_key}, \code{source}, \code{moa_text},
#'   \code{target_symbol}, \code{target_uniprot} and
#'   \code{mechanism_comment} (ChEMBL's free-text note on the mechanism
#'   record, available only when the ChEMBL result was fetched with
#'   \code{fields = "all"}; \code{NA} otherwise).
#' @examples
#' chembl <- data.frame(
#'     QueryIDs = c("CHEMBL941", "CHEMBL941"), chembl_id = "CHEMBL941",
#'     Drug_Name = "IMATINIB", MOA = "Bcr/Abl fusion protein inhibitor",
#'     Action_Type = "INHIBITOR", UniProt_ID = c("P00519", "P11274"))
#' broad <- data.frame(
#'     QueryIDs = "imatinib", pert_iname = "imatinib",
#'     InChIKey = "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
#'     moa = "Bcr-Abl kinase inhibitor", target_gene = "ABL1")
#' ## one MOA, two targets - and Broad contributes no rows at all
#' assembleMoaTargets(list(chembl = chembl, broad = broad))
#' @seealso \code{\link{assembleMoaTable}}, \code{\link{listMoaSources}}
#' @export
assembleMoaTargets <- function(results, resolveGeneSymbol = FALSE,
                               taxId = 9606L) {
    emptyOut <- as.data.frame(stats::setNames(
        replicate(length(.dtiMoaTargetCols), character(0), simplify = FALSE),
        .dtiMoaTargetCols))
    if (length(results) == 0L) return(emptyOut)
    ext <- .dtiMoaExtract(results)
    if (is.null(ext)) return(emptyOut)
    ## Sources that name no target contribute no link rows.
    ext <- ext[!(is.na(ext$target_symbol) & is.na(ext$target_uniprot)), ,
               drop = FALSE]
    if (nrow(ext) == 0L) return(emptyOut)
    out <- unique(ext[, .dtiMoaTargetCols, drop = FALSE])
    rownames(out) <- NULL

    if (isTRUE(resolveGeneSymbol)) {
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

#' Query the MOA-carrying sources and assemble their MOA table in one call
#'
#' Convenience wrapper: \code{\link{queryDrugTargets}} followed by
#' \code{\link{assembleMoaTable}}. Call the two separately when you also
#' want the target link (\code{\link{assembleMoaTargets}}) or each
#' source's full native columns, both of which this wrapper discards.
#'
#' @param queryBy a \code{list(molType, idType, ids)} as accepted by
#'   \code{\link{queryDrugTargets}}.
#' @param sources character vector of sources to query; defaults to the
#'   three that carry MOA terms (see \code{\link{listMoaSources}}).
#' @param ... further arguments passed to \code{\link{queryDrugTargets}},
#'   e.g. \code{brhDbPath}, \code{unichemDbPath}, \code{verbose}.
#' @return The \code{data.frame} described in \code{\link{assembleMoaTable}}.
#' @note Querying compounds by \code{idType = "name"} routes through
#'   \code{\link{queryDrugTargets}}'s compound resolution, which reaches
#'   ChEMBL and Open Targets only via a ChEMBL ID - so those two need
#'   \code{unichemDbPath} (see \code{\link{buildUnichemDb}}) when the
#'   query starts from a name. Without it they are skipped and the table
#'   is assembled from whatever remains.
#' @examples
#' \donttest{
#'   queryMoa(list(molType = "cmp", idType = "chembl_id", ids = "CHEMBL941"),
#'            sources = c("chembl", "opentargets"))
#' }
#' @seealso \code{\link{assembleMoaTable}}, \code{\link{assembleMoaTargets}}
#' @export
queryMoa <- function(queryBy, sources = names(.dtiMoaSourceSpec), ...) {
    assembleMoaTable(queryDrugTargets(queryBy, sources = sources, ...))
}

#' Collapse a MOA table to one row per drug
#'
#' Reshapes \code{\link{assembleMoaTable}}'s output column-wise, with a
#' list-column of the distinct MOA terms each drug carries.
#'
#' @param x a \code{data.frame} as returned by \code{\link{assembleMoaTable}}.
#' @return A \code{data.frame} with one row per \code{query_id} x
#'   \code{drug_key}, columns \code{query_id}, \code{drug_key},
#'   \code{drug_name}, \code{n_sources}, \code{n_moa}, and the
#'   list-columns \code{sources}, \code{moa_text} and \code{action_type}.
#' @examples
#' chembl <- data.frame(
#'     QueryIDs = "CHEMBL941", chembl_id = "CHEMBL941",
#'     Drug_Name = "IMATINIB",
#'     MOA = c("Bcr/Abl fusion protein inhibitor",
#'             "Stem cell growth factor receptor inhibitor"),
#'     Action_Type = "INHIBITOR", UniProt_ID = c("P00519", "P10721"))
#' wide <- moaWide(assembleMoaTable(list(chembl = chembl)))
#' wide[, c("drug_key", "drug_name", "n_sources", "n_moa")]
#' wide$moa_text[[1]]
#' @seealso \code{\link{assembleMoaTable}}
#' @export
moaWide <- function(x) {
    .dtiAssertMoaTable(x)
    if (nrow(x) == 0L) {
        out <- data.frame(query_id = character(0), drug_key = character(0),
                          drug_name = character(0), n_sources = integer(0),
                          n_moa = integer(0), stringsAsFactors = FALSE)
        out$sources <- list(); out$moa_text <- list(); out$action_type <- list()
        return(out)
    }
    grp <- paste(x$query_id, x$drug_key, sep = "\r")
    idx <- split(seq_len(nrow(x)), factor(grp, levels = unique(grp)))
    uniqNoNA <- function(v) unique(v[!is.na(v)])
    out <- data.frame(
        query_id  = vapply(idx, function(i) x$query_id[i][1L], character(1)),
        drug_key  = vapply(idx, function(i) x$drug_key[i][1L], character(1)),
        drug_name = vapply(idx, function(i) {
            nm <- uniqNoNA(x$drug_name[i]); if (length(nm)) nm[1L] else NA_character_
        }, character(1)),
        n_sources = vapply(idx, function(i) length(unique(x$source[i])), integer(1)),
        n_moa     = vapply(idx, function(i) length(uniqNoNA(x$moa_text[i])), integer(1)),
        stringsAsFactors = FALSE)
    out$sources <- unname(lapply(idx, function(i) unique(x$source[i])))
    out$moa_text <- unname(lapply(idx, function(i) uniqNoNA(x$moa_text[i])))
    out$action_type <- unname(lapply(idx, function(i) uniqNoNA(x$action_type[i])))
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
