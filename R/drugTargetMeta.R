## =====================================================================
##  drugTargetMeta.R
##  Cross-source drug-target ANNOTATION query dispatcher for the
##  drugTargetInteractions Bioconductor package: given an identifier of
##  (in principle) any recognised type, resolves it to whatever native ID
##  each requested source needs and dispatches to that source's own
##  bidirectional annotation query function (getChemblDrugTarget(),
##  getDgidbDrugTarget(), getOpenTargetsDrugTarget(), ttdTargetAnnot()) -
##  the original motivating goal of the ID-translation layer (see
##  idTranslation.R / unichemAccess.R, both 2026-07-17), tying it to 6 of
##  the already-ported per-source functions.
##
##  PubChem is deliberately NOT one of these 4: getPubchemDrugTarget()
##  returns raw bioassay measurements, not curated annotations, and so
##  doesn't belong in this annotation-only dispatcher/combiner - see
##  getChemblBioassay()/listBioassayFields() (apiAccess.R) for the
##  bioassay-track counterpart, and the "Bioassay Queries" vignette
##  section.
##
##  Named queryDrugTargets(), not getDrugTarget() - that name is already
##  taken by an older, narrower function in
##  drugTargetAnnotations_Fct.R that queries a pre-generated flat file,
##  a different mechanism entirely.
##
##  Deliberately returns a named list of per-source results, not one
##  unified table: cross-source schema harmonization (e.g. a proper S4
##  container) was explicitly deferred to a later phase, and forcing
##  unification here would be premature - this function's job is ID
##  resolution + dispatch, not schema design.
##
##  Two canonical ID vocabularies, distinct from (but resolved via) the
##  source-native ones each per-source function expects:
##   - gene/protein: "symbol" (HGNC), "uniprot" (accession), "ensembl"
##     (gene ID) - resolved via getUniprotMapping() (idTranslation.R),
##     which already accepts arbitrary UniProt `from`/`to` names, so
##     this vocabulary is easy to extend later without touching the
##     resolution mechanism itself.
##   - compound: "name", "chembl_id", "pubchem_id", "drugbank_id",
##     "chebi_id" - resolved via getUnichemMapping() (unichemAccess.R),
##     equally extensible to any source unichemAccess.R's local database
##     carries. Compound *names* are not a UniChem concept (it is
##     purely structure-based), so a name is either passed straight
##     through (to the 4 sources that accept names natively) or, when a
##     name-only source needs to be reached FROM a structured ID,
##     resolved via UniChem -> ChEMBL ID -> getChemblMolecule()$pref_name.
## =====================================================================


## ---------------------------------------------------------------------
## Gene/protein-side resolution
## ---------------------------------------------------------------------

#' Canonical gene/protein idType -> UniProt REST database name
#'
#' Separate \code{from}/\code{to} vocabularies because UniProt's ID
#' mapping API is not symmetric for the UniProtKB side specifically -
#' live-confirmed 2026-07-17: \code{from="UniProtKB_AC-ID"} is required
#' to map *out of* an accession, but \code{to="UniProtKB-Swiss-Prot"} is
#' required to map *into* one (the canonical reviewed-entry target);
#' using the wrong one for the wrong direction is a 400 Bad Request, not
#' a graceful empty result.
#' @keywords internal
.dtiGeneIdTypeMapFrom <- c(symbol = "Gene_Name", uniprot = "UniProtKB_AC-ID",
                           ensembl = "Ensembl")
#' @rdname dot-dtiGeneIdTypeMapFrom
#' @keywords internal
.dtiGeneIdTypeMapTo <- c(symbol = "Gene_Name", uniprot = "UniProtKB-Swiss-Prot",
                         ensembl = "Ensembl")

#' Translate gene/protein identifiers to a target canonical type
#'
#' Canonical vocabulary: \code{"symbol"} (HGNC gene symbol),
#' \code{"uniprot"} (UniProt accession), \code{"ensembl"} (Ensembl gene
#' ID). A passthrough (no network call) when \code{idType == to}.
#'
#' UniProt's ID mapping API requires one side of any single mapping call
#' to be a UniProtKB database - live-confirmed 2026-07-17:
#' \code{Gene_Name -> Ensembl} directly is a 400 Bad Request, but
#' \code{Gene_Name -> UniProtKB-Swiss-Prot} and
#' \code{UniProtKB_AC-ID -> Ensembl} both work. So whenever neither
#' \code{idType} nor \code{to} is \code{"uniprot"}, this routes through
#' a UniProt accession as an intermediate (two \code{\link{getUniprotMapping}}
#' calls instead of one) automatically - callers never need to know this.
#'
#' If an input ID maps to more than one entry (e.g. multiple reviewed
#' isoforms), the first one returned is used - this function always
#' returns exactly one resolved value per input id.
#'
#' @param ids character vector of source identifiers.
#' @param idType character(1) one of \code{"symbol"}, \code{"uniprot"},
#'   \code{"ensembl"} - the type of \code{ids}.
#' @param to character(1) one of the same three - the type to resolve to.
#' @param taxId integer(1) NCBI taxonomy ID passed to
#'   \code{\link{getUniprotMapping}} (default 9606L = human).
#' @param verbose logical(1); if TRUE, message progress.
#' @return A character vector the same length as \code{ids}, named by
#'   \code{ids}, with resolved identifiers (or \code{NA} for any \code{id}
#'   that did not resolve).
#' @keywords internal
.resolveGeneIds <- function(ids, idType, to, taxId = 9606L, verbose = FALSE) {
    idType <- match.arg(idType, names(.dtiGeneIdTypeMapFrom))
    to <- match.arg(to, names(.dtiGeneIdTypeMapFrom))
    if (identical(idType, to)) return(stats::setNames(ids, ids))

    if (idType != "uniprot" && to != "uniprot") {
        viaUniprot <- .resolveGeneIds(ids, idType = idType, to = "uniprot",
                                      taxId = taxId, verbose = verbose)
        ok <- !is.na(viaUniprot)
        resolved <- stats::setNames(rep(NA_character_, length(ids)), ids)
        if (any(ok)) {
            final <- .resolveGeneIds(unname(viaUniprot[ok]), idType = "uniprot",
                                     to = to, taxId = taxId, verbose = verbose)
            resolved[ok] <- final[viaUniprot[ok]]
        }
        return(resolved)
    }

    if (verbose) message("Resolving ", length(unique(ids)), " gene ID(s): ",
                         idType, " -> ", to)
    m <- getUniprotMapping(ids, from = .dtiGeneIdTypeMapFrom[[idType]],
                           to = .dtiGeneIdTypeMapTo[[to]], taxId = taxId,
                           verbose = verbose)
    resolved <- m$To[match(ids, m$From)]
    stats::setNames(resolved, ids)
}


## ---------------------------------------------------------------------
## Compound-side resolution
## ---------------------------------------------------------------------

#' Canonical compound idType (excluding \code{"name"}) -> UniChem source name
#' @keywords internal
.dtiCompoundIdTypeMap <- c(chembl_id = "chembl", pubchem_id = "pubchem",
                           drugbank_id = "drugbank", chebi_id = "chebi")

#' Translate compound identifiers to a target canonical type
#'
#' Canonical vocabulary: \code{"name"} (compound/drug name - unreliable
#' due to ambiguity, but what most of the 5 drug-target sources accept
#' natively), plus any structured-ID type in
#' \code{\link{.dtiCompoundIdTypeMap}} (\code{"chembl_id"},
#' \code{"pubchem_id"}, \code{"drugbank_id"}, \code{"chebi_id"} -
#' extensible to any source \code{\link{getUnichemMapping}}'s underlying
#' database carries).
#'
#' UniChem is purely structure-based (see \code{unichemAccess.R}) - it
#' has no concept of a compound name - so a \code{"name"} endpoint on
#' either side is handled outside UniChem:
#' \itemize{
#'   \item structured ID -> name: resolved via UniChem to a ChEMBL ID
#'     (unless already one), then \code{\link{getChemblMolecule}}'s
#'     \code{pref_name} - \code{NA} if ChEMBL has no preferred name for
#'     that entry.
#'   \item name -> structured ID: resolved via PubChem's internal name
#'     resolver (\code{.dtiPubchemNameToCid()}, \code{apiAccess.R}) to a
#'     PubChem CID, then via UniChem if a different target type was
#'     requested.
#' }
#' Two structured ID types resolve via one direct UniChem self-join, no
#' name involved at all.
#'
#' @param ids character vector of source identifiers.
#' @param idType character(1) \code{"name"} or a key of
#'   \code{\link{.dtiCompoundIdTypeMap}}.
#' @param to character(1) same vocabulary - the type to resolve to.
#' @param unichemDbPath character(1) path to a UniChem SQLite (see
#'   \code{\link{buildUnichemDb}}); required whenever a structured ID
#'   type is involved on either side and \code{idType != to}.
#' @param verbose logical(1); if TRUE, message progress.
#' @return A character vector the same length as \code{ids}, named by
#'   \code{ids}, with resolved identifiers (or \code{NA} for any \code{id}
#'   that did not resolve).
#' @keywords internal
.resolveCompoundIds <- function(ids, idType, to, unichemDbPath = NULL, verbose = FALSE) {
    validTypes <- c("name", names(.dtiCompoundIdTypeMap))
    idType <- match.arg(idType, validTypes)
    to <- match.arg(to, validTypes)
    if (identical(idType, to)) return(stats::setNames(ids, ids))

    ## structured -> structured: one direct UniChem self-join
    if (idType != "name" && to != "name") {
        if (is.null(unichemDbPath))
            stop("Resolving '", idType, "' to '", to, "' requires unichemDbPath ",
                 "(see buildUnichemDb()).")
        if (verbose) message("Resolving ", length(unique(ids)), " compound ID(s): ",
                             idType, " -> ", to, " via UniChem")
        m <- getUnichemMapping(ids, from = .dtiCompoundIdTypeMap[[idType]],
                               to = .dtiCompoundIdTypeMap[[to]], dbPath = unichemDbPath)
        resolved <- m$To[match(ids, m$From)]
        return(stats::setNames(resolved, ids))
    }

    ## structured -> name: via ChEMBL ID + getChemblMolecule()$pref_name
    if (idType != "name" && to == "name") {
        chemblIds <- if (identical(idType, "chembl_id")) {
            stats::setNames(ids, ids)
        } else {
            .resolveCompoundIds(ids, idType = idType, to = "chembl_id",
                                unichemDbPath = unichemDbPath, verbose = verbose)
        }
        ok <- !is.na(chemblIds)
        resolved <- stats::setNames(rep(NA_character_, length(ids)), ids)
        if (any(ok)) {
            mol <- getChemblMolecule(unique(unname(chemblIds[ok])))
            nameFor <- stats::setNames(mol$pref_name, mol$chembl_id)
            resolved[ok] <- nameFor[chemblIds[ok]]
        }
        return(resolved)
    }

    ## name -> structured: via PubChem's internal name resolver
    if (identical(idType, "name") && to != "name") {
        if (verbose) message("Resolving ", length(unique(ids)), " compound name(s) ",
                             "-> PubChem CID")
        cids <- character(length(ids))
        for (i in seq_along(ids)) {
            cids[i] <- .dtiPubchemNameToCid(ids[i])
            Sys.sleep(0.2)
        }
        cids <- stats::setNames(cids, ids)
        if (identical(to, "pubchem_id")) return(cids)
        ok <- !is.na(cids)
        resolved <- stats::setNames(rep(NA_character_, length(ids)), ids)
        if (any(ok)) {
            final <- .resolveCompoundIds(unname(cids[ok]), idType = "pubchem_id",
                                         to = to, unichemDbPath = unichemDbPath,
                                         verbose = verbose)
            resolved[ok] <- final[cids[ok]]
        }
        return(resolved)
    }
}


## ---------------------------------------------------------------------
## Meta-dispatcher
## ---------------------------------------------------------------------

#' All sources \code{\link{queryDrugTargets}} knows how to dispatch to,
#' and the canonical resolution target each one's native \code{queryBy}
#' interface needs. Gene/protein side all resolve to a symbol - the one
#' type every non-ChEMBL source accepts directly - except ChEMBL, which
#' needs a UniProt accession. Compound side all resolve to a name -
#' efficient in the common case (a name input is a zero-network
#' passthrough - see \code{\link{.resolveCompoundIds}}) - except ChEMBL
#' and Open Targets, which take a ChEMBL ID directly (and, for Open
#' Targets, still work when handed one under \code{idType = "name"},
#' since \code{getOpenTargetsDrugTarget()} passes ChEMBL-shaped strings
#' through unresolved).
#'
#' PubChem is deliberately \strong{not} listed here: its REST accessors
#' (\code{getPubchemDrugs()}/\code{getPubchemTargets()}/
#' \code{getPubchemDrugTarget()}) return raw bioassay measurements, not
#' curated drug-target annotations the way every other source here does
#' - a different data type that doesn't belong in this dispatcher or in
#' \code{\link{combineDrugTargets}}'s column harmonization. See the
#' "Bioassay Queries" vignette section, \code{\link{getChemblBioassay}}
#' and \code{\link{listBioassayFields}} for the bioassay-track
#' counterpart.
#' @keywords internal
.dtiMetaSources <- list(
    chembl      = list(gene = "uniprot",  cmp = "chembl_id"),
    dgidb       = list(gene = "symbol",   cmp = "name"),
    opentargets = list(gene = "symbol",   cmp = "chembl_id"),
    ttd         = list(gene = "symbol",   cmp = "name"),
    broad       = list(gene = "symbol",   cmp = "name"),
    gtopdb      = list(gene = "symbol",   cmp = "name")
)

#' Build the native \code{queryBy} list a given source's function expects
#' @keywords internal
.dtiMetaQueryBy <- function(src, isGene, ids) {
    switch(src,
        chembl = if (isGene) list(molType = "protein", idType = "Uniprot", ids = ids)
                 else list(molType = "cmp", idType = "chembl_id", ids = ids),
        dgidb = if (isGene) list(molType = "gene", idType = "symbol", ids = ids)
                else list(molType = "cmp", idType = "name", ids = ids),
        opentargets = if (isGene) list(molType = "gene", idType = "symbol", ids = ids)
                      else list(molType = "cmp", idType = "name", ids = ids),
        ttd = if (isGene) list(molType = "protein", idType = "symbol", ids = ids)
              else list(molType = "cmp", idType = "name", ids = ids),
        broad = if (isGene) list(molType = "protein", idType = "symbol", ids = ids)
                else list(molType = "cmp", idType = "name", ids = ids),
        gtopdb = if (isGene) list(molType = "protein", idType = "symbol", ids = ids)
                 else list(molType = "cmp", idType = "name", ids = ids))
}

#' Query drug-target interactions across multiple sources from any
#' recognised starting identifier
#'
#' The tying-together piece of the ID-translation layer (see
#' \code{idTranslation.R}, \code{unichemAccess.R}): resolves
#' \code{queryBy$ids} to whatever native identifier each requested
#' source needs (see \code{\link{.resolveGeneIds}}/
#' \code{\link{.resolveCompoundIds}} for how), then dispatches to that
#' source's own bidirectional \emph{annotation} query function
#' (\code{\link{getChemblDrugTarget}}, \code{\link{getDgidbDrugTarget}},
#' \code{\link{getOpenTargetsDrugTarget}}, \code{\link{ttdTargetAnnot}},
#' \code{\link{broadRepurposingHubAnnot}}, \code{\link{gtoPdbTargetAnnot}}).
#' Named \code{queryDrugTargets()}, not \code{getDrugTarget()} - that
#' name is already taken by an older, unrelated function in
#' \code{drugTargetAnnotations_Fct.R}.
#'
#' PubChem is intentionally not one of the dispatchable \code{sources}:
#' its REST accessors return raw bioassay measurements, not curated
#' drug-target annotations, a different data type from what this
#' function and \code{\link{combineDrugTargets}} harmonize. Use
#' \code{\link{getPubchemDrugTarget}} directly, or
#' \code{\link{getChemblBioassay}} for ChEMBL's own bioassay data - see
#' the "Bioassay Queries" vignette section.
#'
#' Returns a named list, one element per successfully-queried source -
#' \strong{not} one unified table. Cross-source schema harmonization (e.g.
#' a proper S4 container) is a separate, deferred phase; this function's
#' job is ID resolution and dispatch, not schema design. A source
#' contributes no list element
#' at all (rather than an empty data.frame) when none of \code{queryBy$ids}
#' resolved to that source's native ID type, or when the source's own
#' function raised an error (e.g. a required local database path was not
#' supplied) - set \code{verbose = TRUE} to see why.
#'
#' @param queryBy named list with components \code{molType}, \code{idType}
#'   and \code{ids}. \code{molType} is \code{"gene"} or \code{"protein"}
#'   (accepted interchangeably) for the target -> drug direction, or
#'   \code{"cmp"} for drug -> target. \code{idType} is one of
#'   \code{"symbol"}, \code{"uniprot"}, \code{"ensembl"} for the gene/
#'   protein side (see \code{\link{.resolveGeneIds}}), or one of
#'   \code{"name"}, \code{"chembl_id"}, \code{"pubchem_id"},
#'   \code{"drugbank_id"}, \code{"chebi_id"} for the compound side (see
#'   \code{\link{.resolveCompoundIds}}) - not each source's own native
#'   vocabulary, which this function translates to internally.
#' @param sources character vector, any of \code{"chembl"},
#'   \code{"dgidb"}, \code{"opentargets"}, \code{"ttd"}, \code{"broad"},
#'   \code{"gtopdb"} (default: all six annotation sources; PubChem is not
#'   included here, see Details).
#' @param ttdDbPath character(1) path to a local TTD SQLite (see
#'   \code{\link{buildTtdDb}}); required if \code{"ttd"} is in
#'   \code{sources}. Not built automatically - a TTD build is a real,
#'   deliberate operation, not something to trigger silently from inside
#'   a dispatcher.
#' @param brhDbPath character(1) path to a local Broad Repurposing Hub
#'   SQLite (see \code{\link{buildBroadRepurposingHubDb}}); required if
#'   \code{"broad"} is in \code{sources}. Not built automatically, same
#'   rationale as \code{ttdDbPath}.
#' @param gtoPdbDbPath character(1) path to a local GtoPdb SQLite (see
#'   \code{\link{buildGtoPdbDb}}); required if \code{"gtopdb"} is in
#'   \code{sources}. Not built automatically, same rationale as
#'   \code{ttdDbPath}.
#' @param unichemDbPath character(1) path to a local UniChem SQLite (see
#'   \code{\link{buildUnichemDb}}); required whenever compound-side
#'   resolution needs it (structured-ID-to-structured-ID or
#'   structured-ID-to-name translation - see
#'   \code{\link{.resolveCompoundIds}}). Not built automatically -
#'   \code{\link{buildUnichemDb}} takes on the order of an hour.
#' @param taxId integer(1) NCBI taxonomy ID for gene/protein-side
#'   resolution (default 9606L = human); passed to
#'   \code{\link{getUniprotMapping}}.
#' @param verbose logical(1); if TRUE, message resolution and per-source
#'   progress, and why a source was skipped.
#' @param ... additional arguments passed through to each dispatched
#'   source function (e.g. \code{expand} for
#'   \code{\link{getOpenTargetsDrugTarget}}).
#' @return A named list, one element (a \code{data.frame}) per source
#'   that returned a result. \code{attr(result, "resolved")} holds the
#'   per-source ID-resolution vectors (named by the original
#'   \code{queryBy$ids}, as returned by \code{\link{.resolveGeneIds}}/
#'   \code{\link{.resolveCompoundIds}}) for tracing which input ID led to
#'   which rows.
#' @examples
#' \donttest{
#'   ## target -> drug, starting from an Ensembl gene ID, all live sources
#'   res <- queryDrugTargets(
#'     list(molType = "gene", idType = "ensembl", ids = "ENSG00000077782"),
#'     sources = c("chembl", "dgidb", "opentargets"))
#'   names(res)
#'   res$chembl
#'
#'   ## drug -> target, starting from a DrugBank ID (needs a UniChem SQLite)
#'   res2 <- queryDrugTargets(
#'     list(molType = "cmp", idType = "drugbank_id", ids = "DB00945"),
#'     sources = "chembl", unichemDbPath = buildUnichemDb(rerun = FALSE))
#' }
#' @seealso \code{\link{getChemblDrugTarget}}, \code{\link{getDgidbDrugTarget}},
#'   \code{\link{getOpenTargetsDrugTarget}}, \code{\link{ttdTargetAnnot}},
#'   \code{\link{broadRepurposingHubAnnot}}, \code{\link{gtoPdbTargetAnnot}},
#'   \code{\link{getPubchemDrugTarget}}, \code{\link{getChemblBioassay}},
#'   \code{\link{getUniprotMapping}}, \code{\link{getUnichemMapping}}
#' @export
queryDrugTargets <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                             sources = names(.dtiMetaSources), ttdDbPath = NULL,
                             brhDbPath = NULL, gtoPdbDbPath = NULL, unichemDbPath = NULL,
                             taxId = 9606L, verbose = FALSE, ...) {
    if (!identical(names(queryBy), c("molType", "idType", "ids"))) {
        stop(
            "All three list components in 'queryBy' (named: 'molType',",
            " 'idType' and 'ids') need to be present."
        )
    }
    if (any(vapply(queryBy, length, integer(1)) == 0)) {
        stop(
            "All components in 'queryBy' list need to be populated with ",
            "corresponding character vectors."
        )
    }
    if (!all(sources %in% names(.dtiMetaSources)))
        stop("'sources' must be one or more of: ",
             paste(names(.dtiMetaSources), collapse = ", "))
    isGene <- identical(queryBy$molType, "gene") || identical(queryBy$molType, "protein")
    if (!isGene && !identical(queryBy$molType, "cmp"))
        stop("queryBy$molType must be \"gene\"/\"protein\" or \"cmp\"")

    out <- list()
    resolvedBySource <- list()
    for (src in sources) {
        needed <- if (isGene) .dtiMetaSources[[src]]$gene else .dtiMetaSources[[src]]$cmp

        result <- tryCatch({
            resolved <- if (isGene) {
                .resolveGeneIds(queryBy$ids, idType = queryBy$idType, to = needed,
                                taxId = taxId, verbose = verbose)
            } else {
                .resolveCompoundIds(queryBy$ids, idType = queryBy$idType, to = needed,
                                    unichemDbPath = unichemDbPath, verbose = verbose)
            }
            resolvedBySource[[src]] <- resolved
            rids <- unique(stats::na.omit(unname(resolved)))
            if (length(rids) == 0L) {
                if (verbose) message("queryDrugTargets: no IDs resolved for '", src, "'")
                return(NULL)
            }
            qb <- .dtiMetaQueryBy(src, isGene, rids)
            if (verbose) message("queryDrugTargets: querying '", src, "' with ",
                                 length(rids), " resolved ID(s)")
            switch(src,
                chembl      = getChemblDrugTarget(qb, verbose = verbose, ...),
                dgidb       = getDgidbDrugTarget(qb, verbose = verbose, ...),
                opentargets = getOpenTargetsDrugTarget(qb, verbose = verbose, ...),
                ttd         = {
                    if (is.null(ttdDbPath))
                        stop("'ttd' requires ttdDbPath (see buildTtdDb()).")
                    ttdTargetAnnot(qb, ttdDbPath)
                },
                broad       = {
                    if (is.null(brhDbPath))
                        stop("'broad' requires brhDbPath (see buildBroadRepurposingHubDb()).")
                    broadRepurposingHubAnnot(qb, brhDbPath)
                },
                gtopdb      = {
                    if (is.null(gtoPdbDbPath))
                        stop("'gtopdb' requires gtoPdbDbPath (see buildGtoPdbDb()).")
                    gtoPdbTargetAnnot(qb, gtoPdbDbPath)
                })
        }, error = function(e) {
            if (verbose) message("queryDrugTargets: source '", src, "' failed: ",
                                 conditionMessage(e))
            NULL
        })
        if (!is.null(result)) out[[src]] <- result
    }

    attr(out, "resolved") <- resolvedBySource
    out
}


## ---------------------------------------------------------------------
## Combining results across sources (row/column append, not
## harmonization - see the file header for why a proper S4 container
## is a deferred, more principled version of this)
## ---------------------------------------------------------------------

#' Display name per source, used as a constant \code{source} column value
#'
#' PubChem is not included - see \code{\link{.dtiMetaSources}} for why
#' it's excluded from this whole annotation-combining layer.
#' @keywords internal
.dtiCombineSourceLabel <- c(chembl = "ChEMBL", dgidb = "DGIdb",
                            opentargets = "OpenTargets", ttd = "TTD",
                            broad = "Broad Repurposing Hub", gtopdb = "GtoPdb")

#' Canonical combined column -> per-source column name.
#'
#' \code{NA} means that source has no equivalent column at all (not that
#' the column exists but is empty) - currently only ChEMBL's
#' \code{gene_symbol} (ChEMBL's REST output is UniProt-accession-keyed,
#' \code{UniProt_ID}, with no gene symbol column). See
#' \code{\link{combineDrugTargets}}'s \code{resolveGeneSymbol} argument
#' for filling this in on request rather than always paying for it.
#' \code{query_id} and \code{source} are handled separately in
#' \code{\link{combineDrugTargets}} (query_id preferably from
#' \code{queryDrugTargets()}'s \code{"resolved"} attribute so it reflects
#' the user's *original* query token rather than each source's own,
#' possibly-translated \code{QueryIDs}; source is always a constant, see
#' \code{\link{.dtiCombineSourceLabel}}, not looked up per row).
#' @keywords internal
.dtiCombineColMap <- list(
    gene_symbol = c(chembl = NA, dgidb = "gene_name",
                    opentargets = "approved_symbol", ttd = "GeneName",
                    broad = "target_gene", gtopdb = "target_gene"),
    drug_name   = c(chembl = "Drug_Name", dgidb = "drug_name",
                    opentargets = "drug_name", ttd = "DrugName",
                    broad = "pert_iname", gtopdb = "ligandName"),
    action      = c(chembl = "Action_Type", dgidb = "interaction_types",
                    opentargets = "action_type", ttd = "MOA",
                    broad = "moa", gtopdb = "action")
)

#' Combine \code{\link{queryDrugTargets}} results into one table
#'
#' Row-binds a small set of canonical columns across whichever sources
#' are present in \code{results}, mapping each source's own column names
#' to a shared vocabulary (see \code{\link{.dtiCombineColMap}}). This is
#' deliberately just column-name alignment and row append - not
#' deduplication, not cross-source identity resolution (the same
#' compound or target appearing under different native IDs in different
#' sources is not merged). A more principled harmonized container (an S4
#' class, analysed but not yet built) is a separate, deferred phase;
#' this function is the "just append them" version to use in the
#' meantime, and every
#' source's full original data remains available unchanged in
#' \code{results} itself.
#'
#' \code{action} is a best-effort common label, not a perfectly aligned
#' concept: ChEMBL's \code{Action_Type} and Open Targets'
#' \code{action_type} are categorical mechanism labels (e.g.
#' \code{"INHIBITOR"}); TTD's \code{MOA}, DGIdb's \code{interaction_types},
#' the Broad Repurposing Hub's \code{moa} (free text, e.g.
#' \code{"fgfr inhibitor"}, not a controlled vocabulary) and GtoPdb's
#' \code{action} (e.g. \code{"Inhibition"}) are similar. PubChem is not
#' one of the
#' combinable sources at all (see \code{\link{.dtiMetaSources}}): its
#' bioactivity data has no mechanism-of-action concept to align to in
#' the first place.
#'
#' @param results a named list as returned by \code{\link{queryDrugTargets}}
#'   (or any similarly-shaped named list of per-source data.frames -
#'   \code{query_id} falls back to each source's own \code{QueryIDs}
#'   column when \code{results} has no \code{"resolved"} attribute).
#' @param columns character vector of canonical columns to include, any
#'   of \code{"query_id"}, \code{"gene_symbol"}, \code{"drug_name"},
#'   \code{"action"}, \code{"source"} (default: all five).
#' @param resolveGeneSymbol logical(1); if \code{TRUE}, fill ChEMBL's
#'   otherwise-\code{NA} \code{gene_symbol} by resolving its
#'   \code{UniProt_ID} column via \code{\link{getUniprotMapping}} (one
#'   extra network round trip - default \code{FALSE} so
#'   \code{combineDrugTargets()} is network-free by default when
#'   \code{results} already has everything it needs). Only has an effect
#'   when \code{"gene_symbol"} is also in \code{columns}.
#' @param taxId integer(1) passed to \code{\link{getUniprotMapping}} when
#'   \code{resolveGeneSymbol = TRUE} (default 9606L = human).
#' @return A single \code{data.frame} with columns \code{columns}, one
#'   row per row of every source in \code{results}.
#' @examples
#' \donttest{
#'   res <- queryDrugTargets(list(molType = "gene", idType = "symbol", ids = "FGFR1"),
#'                           sources = c("chembl", "dgidb", "opentargets"))
#'   combineDrugTargets(res)
#'   combineDrugTargets(res, resolveGeneSymbol = TRUE)  ## fills ChEMBL's gene_symbol too
#' }
#' @seealso \code{\link{queryDrugTargets}}
#' @export
combineDrugTargets <- function(results,
                               columns = c("query_id", "gene_symbol", "drug_name",
                                          "action", "source"),
                               resolveGeneSymbol = FALSE, taxId = 9606L) {
    columns <- match.arg(columns, c("query_id", "gene_symbol", "drug_name",
                                    "action", "source"), several.ok = TRUE)
    if (length(results) == 0L)
        return(as.data.frame(stats::setNames(
            replicate(length(columns), character(0), simplify = FALSE), columns)))

    unsupported <- setdiff(names(results), names(.dtiCombineSourceLabel))
    if (length(unsupported))
        stop("combineDrugTargets() does not recognise source(s): ",
             paste(unsupported, collapse = ", "), ". 'results' must be named ",
             "using the same source keys queryDrugTargets() uses: ",
             paste(names(.dtiCombineSourceLabel), collapse = ", "), ".")

    resolvedAttr <- attr(results, "resolved")
    needCols <- union(columns, if (isTRUE(resolveGeneSymbol)) "gene_symbol" else character(0))

    rows <- lapply(names(results), function(src) {
        df <- results[[src]]
        out <- data.frame(row.names = seq_len(nrow(df)))
        out$.source <- .dtiCombineSourceLabel[[src]]  ## always tracked internally
        if ("query_id" %in% needCols) {
            out$query_id <- if (!is.null(resolvedAttr[[src]])) {
                orig <- names(resolvedAttr[[src]])
                names(orig) <- unname(resolvedAttr[[src]])
                unname(orig[df$QueryIDs])
            } else {
                df$QueryIDs
            }
        }
        for (col in intersect(needCols, names(.dtiCombineColMap))) {
            srcCol <- .dtiCombineColMap[[col]][[src]]
            if (!is.na(srcCol) && !srcCol %in% names(df))
                stop("combineDrugTargets(): expected column '", srcCol, "' not found ",
                     "in results$", src, " - is this really a ", src, " result from ",
                     "queryDrugTargets()?")
            out[[col]] <- if (is.na(srcCol)) NA_character_ else as.character(df[[srcCol]])
        }
        out
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL

    if (isTRUE(resolveGeneSymbol) && "chembl" %in% names(results)) {
        ## rows for a given source were rbound in results$chembl's own row
        ## order, so the ChEMBL block's positions line up 1:1 with
        ## results$chembl$UniProt_ID - no id-matching needed, just a
        ## direct positional assignment (length-safe, unlike an is.na()
        ## mask which could silently misalign if some other column
        ## already happened to be NA for a non-ChEMBL reason).
        chemblIdx <- which(out$.source == "ChEMBL")
        if (length(chemblIdx) > 0L) {
            upIds <- results$chembl$UniProt_ID
            resolved <- .resolveGeneIds(stats::na.omit(unique(upIds)),
                                        idType = "uniprot", to = "symbol", taxId = taxId)
            out$gene_symbol[chemblIdx] <- unname(resolved[upIds])
        }
    }

    if ("source" %in% columns) out$source <- out$.source
    out$.source <- NULL
    out[, columns, drop = FALSE]
}
