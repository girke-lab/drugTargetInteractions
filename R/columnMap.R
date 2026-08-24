## =====================================================================
##  columnMap.R
##  The column mapping table used when per-source tables are appended
##  row-wise into one table.
##
##  Two separate problems arise when stacking the per-source tables:
##
##    1. The same content appears under different titles - a drug name is
##       `Drug_Name` in ChEMBL, `drug_name` in DGIdb and Open Targets,
##       `DrugName` in TTD, `pert_iname` in the Broad Repurposing Hub and
##       `ligandName` in GtoPdb. These have to be aligned onto one column.
##    2. A column exists in only one source - `disease_name` in Open
##       Targets, `Smiles` in TTD. These are appended as their own
##       columns and are NA for every other source's rows.
##
##  Problem 1 is what this file's mapping table solves. It was previously
##  an internal, hardcoded list (.dtiCombineColMap in drugTargetMeta.R)
##  covering three canonical columns and reachable only from
##  combineDrugTargets(); the curated map below is now the single source
##  of truth, and .dtiCombineColMap is derived from it so the two can
##  never drift apart.
##
##  Why the map is curated rather than generated. Matching column names
##  after normalising case and punctuation finds some groups on its own,
##  but it cannot be trusted unsupervised. Checked against a real
##  four-source genome-wide build it proposes
##
##      opentargets:drug_id | ttd:DrugID
##
##  which are not the same identifier space at all - Open Targets'
##  `drug_id` holds a ChEMBL id (CHEMBL1200797), TTD's `DrugID` holds a
##  TTD-internal id (D03DSN) - so aligning them would fabricate identity.
##  It equally misses real groups whose names share nothing, such as
##  `approved_symbol`/`gene_name`/`GeneName`. Name matching therefore
##  produces a draft to review, never a mapping to apply: proposals come
##  back marked inactive and only take effect once a person turns them on.
## =====================================================================


## ---------------------------------------------------------------------
## The curated map
## ---------------------------------------------------------------------
## One row per (canonical column, source). NA in the `column` slot means
## that source has no equivalent column at all - not that the column
## exists and is empty. Only ChEMBL's gene_symbol is like that today
## (ChEMBL's REST output is UniProt-accession-keyed, with no gene symbol
## column); every other absent pair is simply left out of the table.
##
## A native column may legitimately feed more than one canonical column.
## The Broad Repurposing Hub's `moa` is the only such case: it holds free
## text of the form "fgfr inhibitor", which serves both as the coarse
## `action` label it has always been mapped to and as the finer
## `mechanism` descriptor alongside ChEMBL's `MOA`.

#' Curated per-source column mapping, long format
#'
#' The single source of truth for column alignment; see this file's
#' header for why it is curated and not generated.
#' @keywords internal
#' @noRd
.dtiCuratedColumnMap <- local({
    row <- function(canonical, ...) {
        cols <- c(...)
        data.frame(canonical = canonical, source = names(cols),
                   column = unname(cols), stringsAsFactors = FALSE)
    }
    do.call(rbind, list(
        ## Gene symbol as each source reports it. This is what the source
        ## called the gene, which is not always what HGNC calls it today -
        ## a genome-wide table's own `symbol` column is the authoritative
        ## one, and the two are kept side by side so disagreements stay
        ## visible instead of being silently reconciled.
        row("gene_symbol",
            chembl = NA_character_, dgidb = "gene_name",
            opentargets = "approved_symbol", ttd = "GeneName",
            broad = "target_gene", gtopdb = "target_gene"),
        row("drug_name",
            chembl = "Drug_Name", dgidb = "drug_name",
            opentargets = "drug_name", ttd = "DrugName",
            broad = "pert_iname", gtopdb = "ligandName"),
        ## Coarse action label. Vocabularies differ in case and wording
        ## (INHIBITOR / inhibitor / Inhibitor) but the concept is shared.
        row("action",
            chembl = "Action_Type", dgidb = "interaction_types",
            opentargets = "action_type", ttd = "MOA",
            broad = "moa", gtopdb = "action"),
        ## Free-text mechanism naming the target as well as the action
        ## ("26S proteasome inhibitor"), as opposed to `action`'s bare
        ## label. TTD's MOA is a bare label and belongs under `action`,
        ## not here.
        row("mechanism",
            chembl = "MOA", opentargets = "mechanism_of_action",
            broad = "moa"),
        row("target_uniprot",
            chembl = "UniProt_ID", ttd = "Uniprot_acc"),
        row("pubchem_cid",
            chembl = "PubChem_CID", ttd = "PubChem_CID",
            broad = "pubchem_cid"),
        ## Furthest stage of development reached. The concept is shared;
        ## the vocabularies are not interchangeable - ChEMBL reports a
        ## number (-1 to 4), Open Targets a token (PHASE_3, APPROVAL),
        ## TTD and the Broad Hub a phrase (Approved, Phase 1). Read this
        ## column together with `source`.
        row("max_phase",
            chembl = "Max_Phase", opentargets = "max_clinical_stage",
            ttd = "Highest_status", broad = "clinical_phase"),
        ## Diseases the drug is indicated for, as a "; "-joined list.
        ## Each source decorates the names differently: ChEMBL prefixes
        ## MeSH ids, TTD appends ICD-11 codes and a per-disease status,
        ## Open Targets gives bare names.
        row("indication",
            chembl = "Mesh_Indication", opentargets = "disease_name",
            ttd = "Indication", broad = "indication")
    ))
})

#' Columns the map deliberately leaves alone
#'
#' \code{hgnc_id} and \code{compound_chembl_id} are cross-source keys
#' resolved by \code{\link{addCommonIds}} against HGNC and UniChem, with
#' ambiguity rules a column mapping cannot express (116 UniProt
#' accessions name more than one gene; DGIdb's compound ids carry a
#' \code{chembl:}/\code{ncit:}/\code{rxcui:} prefix that has to be read
#' before the value means anything). Mapping a native column onto them
#' would produce a worse answer than asking for the key, so the map
#' stays out of their way.
#' @keywords internal
#' @noRd
.dtiKeyOwnedCols <- c("hgnc_id", "compound_chembl_id")

#' Identity columns a genome-wide build carries, emitted verbatim
#'
#' \code{\link{buildGenomeWideDrugTargetTable}} tags every row with the
#' HGNC gene it queried from, so these are correct by construction and
#' are never re-derived. \code{QueryIDs} is the identifier actually sent
#' to that source.
#' @keywords internal
#' @noRd
.dtiIdentityCols <- c("hgnc_id", "symbol", "ensembl_gene_id", "QueryIDs")

#' The map columns, in output order
#' @keywords internal
#' @noRd
.dtiColumnMapCols <- c("canonical", "source", "column", "status", "active")


## ---------------------------------------------------------------------
## Public accessor
## ---------------------------------------------------------------------

#' Column mapping used to append per-source tables into one table
#'
#' Every drug-target source names its columns differently, so stacking
#' their tables row-wise requires knowing which of them hold the same
#' content. This function returns that mapping as an ordinary
#' data.frame, one row per (canonical column, source), which you can
#' inspect, edit, write to a file and pass back to
#' \code{\link{combineGenomeWideDrugTargets}}.
#'
#' Called with no arguments it returns the curated mapping the package
#' ships. Given a result list it additionally inspects the columns those
#' tables actually carry and reports two further kinds of row:
#' \code{"proposed"} groups it found by matching column names after
#' normalising case and punctuation, and \code{"unmatched"} columns that
#' belong to no group at all. Proposed rows come back with
#' \code{active = FALSE} and have no effect until you set them
#' \code{TRUE} - column names that look alike are not always the same
#' content, so nothing detected this way is applied on your behalf.
#'
#' \code{hgnc_id} and \code{compound_chembl_id} are not part of the
#' mapping. They are cross-source keys that \code{\link{addCommonIds}}
#' resolves against HGNC and UniChem, using ambiguity rules a column
#' mapping cannot express.
#'
#' @param results optional named list of per-source data.frames, as
#'   returned by \code{\link{buildGenomeWideDrugTargetTable}} or
#'   \code{\link{queryDrugTargets}}. When supplied, proposed and
#'   unmatched rows are added for the columns those tables carry.
#' @param propose logical(1); set \code{FALSE} to report only the
#'   unmatched columns of \code{results} and skip name-based proposals.
#'   Ignored when \code{results} is \code{NULL}.
#' @return A \code{data.frame} with columns \code{canonical} (the shared
#'   column name; \code{NA} for unmatched rows), \code{source},
#'   \code{column} (that source's own column name), \code{status}
#'   (\code{"curated"}, \code{"proposed"} or \code{"unmatched"}) and
#'   \code{active} (whether \code{\link{combineGenomeWideDrugTargets}}
#'   will use the row).
#' @examples
#' \donttest{
#'   drugTargetColumnMap()
#'
#'   ## Review what a real build carries beyond the curated map
#'   res <- buildGenomeWideDrugTargetTable(
#'     hgncTable = getHgncGeneTable()[1:3, ],
#'     sources = c("dgidb", "opentargets"), outDir = tempfile("dti_"))
#'   m <- drugTargetColumnMap(res)
#'   subset(m, status != "curated")
#'
#'   ## Accept one proposal and align a group of your own
#'   m$active[m$canonical %in% "druglikeness"] <- TRUE
#'   combineGenomeWideDrugTargets(res, colMap = m)
#' }
#' @seealso \code{\link{combineGenomeWideDrugTargets}},
#'   \code{\link{addCommonIds}}
#' @export
drugTargetColumnMap <- function(results = NULL, propose = TRUE) {
    curated <- .dtiCuratedColumnMap
    curated$status <- "curated"
    curated$active <- TRUE
    curated <- curated[, .dtiColumnMapCols, drop = FALSE]
    if (is.null(results)) {
        rownames(curated) <- NULL
        return(curated)
    }
    if (!is.list(results) || is.data.frame(results))
        stop("'results' must be a named list of per-source data.frames, as ",
             "returned by buildGenomeWideDrugTargetTable() or ",
             "queryDrugTargets().")

    ## Every (source, column) the tables actually carry, minus the ones
    ## no mapping should touch: the identity block a genome-wide build
    ## adds, and the keys addCommonIds() owns.
    present <- do.call(rbind, c(list(NULL), lapply(names(results), function(src) {
        df <- results[[src]]
        if (!is.data.frame(df)) return(NULL)
        keep <- setdiff(names(df), c(.dtiIdentityCols, .dtiKeyOwnedCols))
        if (length(keep) == 0L) return(NULL)
        data.frame(source = src, column = keep, stringsAsFactors = FALSE)
    })))
    if (is.null(present)) {
        rownames(curated) <- NULL
        return(curated)
    }

    ## A column is already accounted for if the curated map names it for
    ## that same source.
    curatedPair <- paste(curated$source, curated$column, sep = "\r")
    present <- present[!paste(present$source, present$column, sep = "\r") %in%
                       curatedPair, , drop = FALSE]

    proposed <- NULL
    if (isTRUE(propose) && nrow(present) > 0L) {
        norm <- tolower(gsub("[^A-Za-z0-9]", "", present$column))
        ## A group needs at least two sources to be worth aligning; a name
        ## appearing twice within one source is not a cross-source group.
        nSrc <- vapply(split(present$source, norm),
                       function(s) length(unique(s)), integer(1))
        groups <- names(nSrc)[nSrc > 1L]
        if (length(groups)) {
            take <- norm %in% groups
            proposed <- data.frame(
                canonical = norm[take], source = present$source[take],
                column = present$column[take], status = "proposed",
                active = FALSE, stringsAsFactors = FALSE)
            present <- present[!take, , drop = FALSE]
        }
    }

    unmatched <- if (nrow(present) > 0L)
        data.frame(canonical = NA_character_, source = present$source,
                   column = present$column, status = "unmatched",
                   active = FALSE, stringsAsFactors = FALSE) else NULL

    out <- do.call(rbind, Filter(Negate(is.null),
                                 list(curated, proposed, unmatched)))
    rownames(out) <- NULL
    out
}


## ---------------------------------------------------------------------
## Validation
## ---------------------------------------------------------------------

#' Check a (possibly user-edited) column map before it is applied
#'
#' Returns the active, usable subset. Errors on the one edit that cannot
#' be resolved sensibly: two different native columns of the same source
#' both claiming one canonical column, which leaves no way to decide
#' which value that column should hold.
#' @keywords internal
#' @noRd
.dtiValidateColumnMap <- function(colMap) {
    if (!is.data.frame(colMap))
        stop("'colMap' must be a data.frame, as returned by ",
             "drugTargetColumnMap().")
    missingCols <- setdiff(c("canonical", "source", "column"), names(colMap))
    if (length(missingCols))
        stop("'colMap' is missing required column(s): ",
             paste(missingCols, collapse = ", "), ". See drugTargetColumnMap().")
    if (is.null(colMap$active)) colMap$active <- TRUE

    keep <- !is.na(colMap$canonical) & .dtiIsTrue(colMap$active) &
            !is.na(colMap$column)
    colMap <- colMap[keep, , drop = FALSE]
    if (nrow(colMap) == 0L) return(colMap)

    clash <- colMap$canonical %in% .dtiKeyOwnedCols
    if (any(clash))
        stop("'colMap' maps column(s) onto ",
             paste(unique(colMap$canonical[clash]), collapse = " and "),
             ", which is resolved by addCommonIds() rather than by column ",
             "mapping - drop those row(s), or rename the canonical column.")

    pair <- paste(colMap$canonical, colMap$source, sep = "\r")
    dup <- unique(pair[duplicated(pair)])
    if (length(dup)) {
        bad <- colMap[pair %in% dup, , drop = FALSE]
        stop("'colMap' is ambiguous: ",
             paste(sprintf("%s <- %s:{%s}",
                           bad$canonical[!duplicated(pair[pair %in% dup])],
                           bad$source[!duplicated(pair[pair %in% dup])],
                           vapply(dup, function(d)
                               paste(colMap$column[pair == d], collapse = ", "),
                               character(1))),
                   collapse = "; "),
             ". One source cannot fill one canonical column from two ",
             "different columns - keep one row per (canonical, source) pair.")
    }
    colMap
}

#' Vectorised isTRUE(), for the map's `active` flag
#' @keywords internal
#' @noRd
.dtiIsTrue <- function(x) !is.na(x) & as.logical(x)
