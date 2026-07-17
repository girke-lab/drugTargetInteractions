## =====================================================================
##  apiAccess.R
##  Live-API drug/target resource access for the drugTargetInteractions
##  Bioconductor package.
##
##  Adds ChEMBL REST query layers alongside the package's existing
##  local-SQLite ChEMBL path. All functions here reach public web APIs
##  at call time (no API key required) and therefore MUST be
##  network-guarded: examples, tests and vignette chunks that call them
##  should be wrapped so that an offline build machine degrades
##  gracefully rather than erroring (see `.dtiApiGET` / `.dtiHasInternet`
##  and the roxygen \donttest{} tags below).
##
##  Naming follows the existing package API (camelCase, `get*` accessors).
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure (not exported)
## ---------------------------------------------------------------------

#' Endpoint registry
#'
#' Central place for the base URLs so a mirror or a pinned proxy can be
#' swapped in one location. Not exported.
#' @keywords internal
.dtiEndpoints <- function() {
    list(
        chembl      = "https://www.ebi.ac.uk/chembl/api/data",
        opentargets = "https://api.platform.opentargets.org/api/v4/graphql",
        dgidb       = "https://dgidb.org/api/graphql",
        pubchem     = "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        eutils      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    )
}

#' Lightweight connectivity probe
#'
#' Used to guard examples/tests. Returns TRUE only if the host resolves
#' and answers within `timeout` seconds. Never throws.
#' @param url character(1) URL to probe (default Open Targets API root).
#' @param timeout numeric(1) seconds.
#' @return logical(1)
#' @keywords internal
.dtiHasInternet <- function(url = "https://www.ebi.ac.uk", timeout = 5) {
    ok <- tryCatch({
        req <- httr2::request(url)
        req <- httr2::req_timeout(req, timeout)
        req <- httr2::req_method(req, "HEAD")
        resp <- httr2::req_perform(req)
        httr2::resp_status(resp) < 500
    }, error = function(e) FALSE)
    isTRUE(ok)
}

#' Perform a GET against a REST endpoint with retry + polite throttling
#'
#' Wraps httr2 with a client-side rate limit (a token-bucket throttle,
#' scoped per-host so different APIs don't share a budget), a bounded
#' exponential back-off on transient failures, an explicit timeout, a
#' descriptive User-Agent, and JSON parsing. On any failure it returns
#' NULL (so callers can degrade gracefully) unless `hardStop = TRUE`.
#' The throttle is a courtesy cap, not a guarantee: none of these APIs
#' publish a hard rate limit, so 429s are still handled gracefully via
#' the retry back-off above.
#'
#' @param url character(1) fully-qualified URL.
#' @param query named list of query parameters (optional).
#' @param timeout numeric(1) per-request timeout in seconds.
#' @param maxTries integer(1) total attempts including the first.
#' @param hardStop logical(1) if TRUE, rethrow the error instead of NULL.
#' @return parsed JSON (list) or NULL on failure.
#' @keywords internal
.dtiApiGET <- function(url, query = NULL, timeout = 60L, maxTries = 3L,
                       hardStop = FALSE) {
    out <- tryCatch({
        req <- httr2::request(url)
        req <- httr2::req_headers(req, Accept = "application/json")
        req <- httr2::req_user_agent(
            req, "drugTargetInteractions R package (Bioconductor)")
        req <- httr2::req_timeout(req, timeout)
        if (!is.null(query))
            req <- httr2::req_url_query(req, !!!query)
        ## Client-side courtesy throttle: at most 5 requests/second per
        ## host (realm defaults to the request's hostname).
        req <- httr2::req_throttle(req, rate = 5, fill_time_s = 1)
        ## Retry on transient conditions; respect Retry-After if present.
        req <- httr2::req_retry(
            req, max_tries = maxTries,
            is_transient = function(resp)
                httr2::resp_status(resp) %in% c(429L, 500L, 502L, 503L, 504L))
        resp <- httr2::req_perform(req)
        httr2::resp_body_json(resp, simplifyVector = FALSE)
    }, error = function(e) {
        if (hardStop) stop(e)
        warning("drugTargetInteractions API GET failed for '", url, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
    })
    out
}

#' Split a vector into chunks of at most `size` elements
#' @keywords internal
.dtiChunk <- function(x, size) {
    if (length(x) == 0L) return(list())
    split(x, ceiling(seq_along(x) / size))
}

#' Batched + paginated GET against a ChEMBL-style list resource
#'
#' Many ChEMBL REST list resources (`target`, `mechanism`,
#' `drug_indication`, `molecule`, ...) accept a `<field>__in=id1,id2,...`
#' filter. A single request's URL is capped by the server at a request-
#' line length of ~4KB, which in practice limits a single `__in` filter
#' to roughly 250-400 IDs depending on ID length (`httr2` percent-encodes
#' the separating commas, which eats into that budget) — so for larger ID
#' sets this chunks `ids` into batches of `chunkSize`, and paginates each
#' batch via `offset`/`page_meta$total_count` until it is exhausted.
#'
#' @param url character(1) resource URL, e.g. `.../mechanism.json`.
#' @param filterField character(1) filter query param, e.g.
#'   `"target_chembl_id__in"`.
#' @param ids character vector of IDs to filter on (deduplicated
#'   internally; order does not matter, this only fetches records).
#' @param resultsField character(1) name of the list field in the JSON
#'   response holding the records, e.g. `"mechanisms"`.
#' @param chunkSize integer(1) max IDs per `__in` filter / request.
#' @param extraQuery named list of additional fixed query parameters.
#' @param verbose logical(1); if TRUE, message progress per batch.
#' @return a flat list of record lists (unparsed JSON records), pooled
#'   across all chunks/pages.
#' @keywords internal
.dtiBatchGET <- function(url, filterField, ids, resultsField, chunkSize = 200L,
                         extraQuery = list(), verbose = FALSE) {
    ids <- unique(stats::na.omit(ids))
    if (length(ids) == 0L) return(list())
    allRecs <- list()
    for (chunk in .dtiChunk(ids, chunkSize)) {
        offset <- 0L
        repeat {
            q <- extraQuery
            q[[filterField]] <- paste(chunk, collapse = ",")
            q$limit <- 1000L
            q$offset <- offset
            if (verbose)
                message("ChEMBL batch GET ", basename(url), " ", filterField,
                        " n=", length(chunk), " offset=", offset)
            res <- .dtiApiGET(url, query = q)
            recs <- res[[resultsField]] %||% list()
            if (length(recs) == 0L) break
            allRecs <- c(allRecs, recs)
            offset <- offset + length(recs)
            total <- res$page_meta$total_count %||% offset
            if (offset >= total) break
        }
    }
    allRecs
}

## NULL-coalescing helper. No roxygen doc block: a leading `%` in \name
## trips up Rd's checkRd (\name should not contain !, | or @).
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a


## ---------------------------------------------------------------------
## ChEMBL REST access
## ---------------------------------------------------------------------

#' Retrieve molecule records from the ChEMBL REST API
#'
#' Queries the public ChEMBL REST endpoint for one or more molecules by
#' ChEMBL ID and returns a tidy \code{data.frame} of the most commonly
#' used identifiers, structures and computed properties. Complements the
#' package's local-SQLite path (\code{\link{downloadChemblDb}}) for cases
#' where only a handful of compounds are needed or a local database is not
#' available.
#'
#' @param chemblIds character vector of ChEMBL molecule IDs
#'   (e.g. \code{"CHEMBL25"}).
#' @param verbose logical(1); if TRUE, message progress per batch.
#' @param chunkSize integer(1) max IDs looked up per HTTP request
#'   (see \code{\link{getChemblDrugTarget}} for why this is capped around
#'   200-400 rather than sent in one request). Large \code{chemblIds}
#'   vectors are chunked and paginated automatically.
#' @return A \code{data.frame} with one row per input ID (duplicates and
#'   input order preserved) and columns: \code{chembl_id}, \code{pref_name},
#'   \code{molecule_type}, \code{max_phase}, \code{first_approval},
#'   \code{canonical_smiles}, \code{standard_inchi_key}, \code{mw_freebase},
#'   \code{alogp}, \code{hba}, \code{hbd}, \code{psa}, \code{rtb},
#'   \code{qed_weighted}. IDs that fail to resolve yield a row of \code{NA}s.
#' @examples
#' \donttest{
#'   ## Requires internet access to www.ebi.ac.uk
#'   df <- getChemblMolecule(c("CHEMBL25", "CHEMBL1201585"))
#'   df[, c("chembl_id", "pref_name", "max_phase")]
#' }
#' @seealso \code{\link{getChemblDrugTarget}}, \code{\link{downloadChemblDb}}
#' @importFrom httr2 request req_headers req_user_agent req_timeout req_url_query req_retry req_throttle req_perform resp_body_json resp_status req_method req_body_json
#' @export
getChemblMolecule <- function(chemblIds, verbose = FALSE, chunkSize = 200L) {
    stopifnot(is.character(chemblIds), length(chemblIds) >= 1L)
    base <- .dtiEndpoints()$chembl
    flat <- function(rec) {
        if (is.null(rec)) return(NULL)
        mp  <- rec$molecule_properties %||% list()
        ms  <- rec$molecule_structures %||% list()
        data.frame(
            chembl_id          = rec$molecule_chembl_id %||% NA_character_,
            pref_name          = rec$pref_name          %||% NA_character_,
            molecule_type      = rec$molecule_type      %||% NA_character_,
            max_phase          = as.numeric(rec$max_phase %||% NA),
            first_approval     = as.integer(rec$first_approval %||% NA),
            canonical_smiles   = ms$canonical_smiles    %||% NA_character_,
            standard_inchi_key = ms$standard_inchi_key  %||% NA_character_,
            mw_freebase        = as.numeric(mp$full_mwt %||% NA),
            alogp              = as.numeric(mp$alogp    %||% NA),
            hba                = as.numeric(mp$hba      %||% NA),
            hbd                = as.numeric(mp$hbd      %||% NA),
            psa                = as.numeric(mp$psa      %||% NA),
            rtb                = as.numeric(mp$rtb      %||% NA),
            qed_weighted       = as.numeric(mp$qed_weighted %||% NA),
            stringsAsFactors   = FALSE)
    }
    recs <- .dtiBatchGET(paste0(base, "/molecule.json"), "molecule_chembl_id__in",
                         chemblIds, "molecules", chunkSize = chunkSize,
                         verbose = verbose)
    found <- lapply(recs, flat)
    foundDF <- if (length(found)) do.call(rbind, found) else flat(list(molecule_chembl_id = "x"))[0, ]
    rownames(foundDF) <- foundDF$chembl_id
    rows <- lapply(chemblIds, function(id) {
        if (id %in% rownames(foundDF)) foundDF[id, ] else flat(list(molecule_chembl_id = id))
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' Search ChEMBL for a target by gene symbol or name
#'
#' Uses the ChEMBL REST \code{target} resource to resolve a query string
#' (gene symbol, protein or preferred name) to ChEMBL target IDs.
#'
#' @param query character(1) search term, e.g. \code{"FGFR1"}.
#' @param organism character(1) optional NCBI organism name filter,
#'   e.g. \code{"Homo sapiens"}; \code{NA} to disable.
#' @param limit integer(1) maximum number of hits to return.
#' @return A \code{data.frame} with columns \code{target_chembl_id},
#'   \code{pref_name}, \code{target_type}, \code{organism}. Empty
#'   \code{data.frame} if nothing matches or the API is unreachable.
#' @examples
#' \donttest{
#'   getChemblTarget("FGFR1", organism = "Homo sapiens")
#' }
#' @seealso \code{\link{getChemblDrugTarget}}
#' @export
getChemblTarget <- function(query, organism = "Homo sapiens", limit = 25L) {
    stopifnot(is.character(query), length(query) == 1L)
    base <- .dtiEndpoints()$chembl
    url  <- paste0(base, "/target/search.json")
    q <- list(q = query, limit = limit)
    res <- .dtiApiGET(url, query = q)
    empty <- data.frame(target_chembl_id = character(0),
                         pref_name = character(0),
                         target_type = character(0),
                         organism = character(0),
                         stringsAsFactors = FALSE)
    if (is.null(res) || is.null(res$targets) || length(res$targets) == 0L)
        return(empty)
    df <- do.call(rbind, lapply(res$targets, function(t) data.frame(
        target_chembl_id = t$target_chembl_id %||% NA_character_,
        pref_name        = t$pref_name        %||% NA_character_,
        target_type      = t$target_type      %||% NA_character_,
        organism         = t$organism         %||% NA_character_,
        stringsAsFactors = FALSE)))
    if (!is.na(organism))
        df <- df[!is.na(df$organism) & df$organism == organism, , drop = FALSE]
    rownames(df) <- NULL
    df
}

#' Retrieve bioactivity measurements from the ChEMBL REST API
#'
#' Pulls activity rows for a given ChEMBL target ID, optionally filtered
#' to a standard activity type (e.g. IC50, Ki). Paginates through the REST
#' \code{activity} resource up to \code{maxRows}.
#'
#' @param targetChemblId character(1), e.g. \code{"CHEMBL1862"}.
#' @param standardType character(1) or NA; e.g. \code{"IC50"}. NA = all.
#' @param maxRows integer(1) cap on returned rows.
#' @param pageSize integer(1) rows per REST page (max 1000).
#' @param verbose logical(1) progress messages.
#' @return A \code{data.frame} with columns \code{molecule_chembl_id},
#'   \code{target_chembl_id}, \code{standard_type}, \code{standard_relation},
#'   \code{standard_value}, \code{standard_units}, \code{pchembl_value},
#'   \code{assay_chembl_id}, \code{assay_description}. Empty on failure.
#' @examples
#' \donttest{
#'   act <- getChemblBioactivities("CHEMBL1862", standardType = "IC50",
#'                                 maxRows = 50)
#'   head(act)
#' }
#' @seealso \code{\link{getChemblTarget}}
#' @export
getChemblBioactivities <- function(targetChemblId, standardType = NA,
                                   maxRows = 1000L, pageSize = 1000L,
                                   verbose = FALSE) {
    stopifnot(is.character(targetChemblId), length(targetChemblId) == 1L)
    base <- .dtiEndpoints()$chembl
    cols <- c("molecule_chembl_id", "target_chembl_id", "standard_type",
              "standard_relation", "standard_value", "standard_units",
              "pchembl_value", "assay_chembl_id", "assay_description")
    empty <- as.data.frame(stats::setNames(
        replicate(length(cols), character(0), simplify = FALSE), cols),
        stringsAsFactors = FALSE)
    q <- list(target_chembl_id = targetChemblId,
              limit = min(pageSize, 1000L), offset = 0L)
    if (!is.na(standardType)) q$standard_type <- standardType
    acc <- list(); got <- 0L; offset <- 0L
    repeat {
        q$offset <- offset
        url <- paste0(base, "/activity.json")
        res <- .dtiApiGET(url, query = q)
        if (is.null(res) || is.null(res$activities) ||
            length(res$activities) == 0L) break
        if (verbose) message("ChEMBL activities offset ", offset,
                             ": +", length(res$activities))
        acc <- c(acc, res$activities)
        got <- got + length(res$activities)
        offset <- offset + length(res$activities)
        total <- res$page_meta$total_count %||% got
        if (got >= maxRows || offset >= total) break
    }
    if (length(acc) == 0L) return(empty)
    acc <- acc[seq_len(min(length(acc), maxRows))]
    df <- do.call(rbind, lapply(acc, function(a) data.frame(
        molecule_chembl_id = a$molecule_chembl_id %||% NA_character_,
        target_chembl_id   = a$target_chembl_id   %||% NA_character_,
        standard_type      = a$standard_type      %||% NA_character_,
        standard_relation  = a$standard_relation  %||% NA_character_,
        standard_value     = as.numeric(a$standard_value %||% NA),
        standard_units     = a$standard_units     %||% NA_character_,
        pchembl_value      = as.numeric(a$pchembl_value %||% NA),
        assay_chembl_id    = a$assay_chembl_id    %||% NA_character_,
        assay_description  = a$assay_description   %||% NA_character_,
        stringsAsFactors   = FALSE)))
    rownames(df) <- NULL
    df
}

#' Resolve the latest ChEMBL release number
#'
#' Queries the ChEMBL REST \code{status} resource for the release currently
#' served (e.g. \code{"ChEMBL_37"}) and returns the bare version number.
#' Used as the default \code{version} for \code{\link{downloadChemblDb}} so
#' the package always pins to the current release unless a specific older
#' version is requested explicitly.
#' @return integer(1) ChEMBL release number.
#' @keywords internal
.chemblLatestVersion <- function() {
    res <- .dtiApiGET(paste0(.dtiEndpoints()$chembl, "/status.json"),
                      hardStop = TRUE)
    v <- res$chembl_db_version %||%
        stop("Could not resolve the latest ChEMBL release from ",
             "the ChEMBL REST status endpoint.")
    as.integer(sub("^ChEMBL_", "", v))
}

#' Query known drug-target annotations via the ChEMBL REST API
#'
#' REST-API equivalent of \code{\link{drugTargetAnnot}} for the package's
#' local ChEMBL SQLite path: given one or more UniProt accessions, returns
#' the drugs annotated (via ChEMBL's \code{drug_mechanism} table) as acting
#' on that target (target -> drug); given one or more ChEMBL molecule IDs,
#' returns the targets annotated for that compound's mechanism of action
#' (drug -> target). Reproduces \code{drugTargetAnnot()}'s join across
#' \code{drug_mechanism}, \code{molecule_dictionary},
#' \code{target_dictionary}/\code{target_components}/
#' \code{component_sequences} and \code{drug_indication} using the ChEMBL
#' REST \code{target}, \code{mechanism}, \code{molecule} and
#' \code{drug_indication} resources, so a local ChEMBL SQLite download is
#' not required.
#'
#' Only ChEMBL's own native identifiers are supported directly:
#' UniProt accession for the target direction, ChEMBL molecule ID for the
#' drug direction. Translating a different starting identifier (Ensembl
#' gene ID, PubChem CID, DrugBank ID, ...) into one of these is left to a
#' separate, source-agnostic ID-translation step (see
#' \code{\link{cmpIdMapping}} for the compound side and
#' \code{\link{getUniprotIDs}}/\code{\link{getParalogs}} for the protein
#' side) so that this function's interface stays uniform across sources.
#'
#' @param queryBy named list with components \code{molType}, \code{idType}
#'   and \code{ids}, matching \code{\link{drugTargetAnnot}}'s interface.
#'   \code{molType = "protein"} with \code{idType} containing
#'   \code{"Uniprot"} (e.g. \code{"Uniprot"}) queries target -> drug by
#'   UniProt accession(s) in \code{ids}. \code{molType = "cmp"} with
#'   \code{idType = "chembl_id"} queries drug -> target by ChEMBL molecule
#'   ID(s) in \code{ids}. Other \code{idType} values (\code{molregno},
#'   \code{PubChem_ID}, \code{DrugBank_ID}) are not yet supported over
#'   REST. \code{ids} can be arbitrarily large: internally the query is
#'   chunked (see \code{chunkSize}) and each chunk paginated, so a query of
#'   e.g. 2000 UniProt accessions or ChEMBL IDs works the same as one.
#' @param verbose logical(1); if TRUE, message progress per batch.
#' @param chunkSize integer(1) max IDs sent per HTTP request via a
#'   \code{<field>__in=id1,id2,...} filter. ChEMBL's REST server caps the
#'   request-line length at ~4KB, which limits a single filter to roughly
#'   250-400 IDs depending on ID length (\code{httr2} percent-encodes the
#'   separating commas, eating into that budget) — 200 leaves comfortable
#'   margin for both short legacy IDs (e.g. \code{"CHEMBL25"}) and longer
#'   modern ones (e.g. \code{"CHEMBL5291677"}). Requests are also
#'   client-side throttled (see \code{\link{.dtiApiGET}}) to be a polite
#'   API citizen; ChEMBL does not publish a hard rate limit.
#' @return A \code{data.frame} with columns \code{QueryIDs},
#'   \code{chembl_id}, \code{Drug_Name}, \code{MOA}, \code{Action_Type},
#'   \code{Max_Phase}, \code{First_Approval}, \code{ChEMBL_TID},
#'   \code{UniProt_ID}, \code{Desc}, \code{Organism},
#'   \code{Mesh_Indication}. Query IDs that return no rows still appear as
#'   a single row with all other fields \code{NA}, so callers can always
#'   confirm which of their input IDs were resolved.
#' @examples
#' \donttest{
#'   ## target -> drug: FGFR1
#'   getChemblDrugTarget(list(molType = "protein", idType = "Uniprot",
#'                            ids = "P11362"))
#'   ## drug -> target: dasatinib
#'   getChemblDrugTarget(list(molType = "cmp", idType = "chembl_id",
#'                            ids = "CHEMBL1421"))
#' }
#' @seealso \code{\link{drugTargetAnnot}}, \code{\link{getChemblMolecule}}
#' @export
getChemblDrugTarget <- function(queryBy = list(molType = NULL, idType = NULL,
                                               ids = NULL),
                                verbose = FALSE, chunkSize = 200L) {
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

    base <- .dtiEndpoints()$chembl
    emptyCols <- c("QueryIDs", "chembl_id", "Drug_Name", "MOA", "Action_Type",
                  "Max_Phase", "First_Approval", "ChEMBL_TID", "UniProt_ID",
                  "Desc", "Organism", "Mesh_Indication")
    emptyDF <- as.data.frame(stats::setNames(
        replicate(length(emptyCols), character(0), simplify = FALSE),
        emptyCols), stringsAsFactors = FALSE)

    isTarget <- identical(queryBy$molType, "protein") &&
        grepl("Uniprot", queryBy$idType[[1]], ignore.case = TRUE)
    isCmp <- identical(queryBy$molType, "cmp") &&
        identical(queryBy$idType, "chembl_id")
    if (!isTarget && !isCmp) {
        stop(
            "getChemblDrugTarget() currently supports only ",
            "queryBy=list(molType=\"protein\", idType=\"Uniprot\", ids=...) ",
            "or queryBy=list(molType=\"cmp\", idType=\"chembl_id\", ids=...). ",
            "Other identifier types (molregno, PubChem_ID, DrugBank_ID) ",
            "require translating to a ChEMBL ID first (e.g. via ",
            "cmpIdMapping()) or the local drugTargetAnnot()."
        )
    }

    ## Resolve UniProt accession -> ChEMBL_TID/Organism/Desc (target -> drug
    ## entry point). One accession can map to several ChEMBL targets
    ## (single-protein plus any protein-family/complex targets it belongs
    ## to), mirroring drugTargetAnnot()'s unfiltered target_components join.
    ## Batched: a returned (possibly multi-component) target is exploded
    ## against every query accession it actually contains, so batching
    ## several accessions into one request reproduces the same per-
    ## accession row set as querying them one at a time.
    .resolveTargetsByAccession <- function(accessions) {
        accessions <- unique(accessions)
        recs <- .dtiBatchGET(paste0(base, "/target.json"),
                             "target_components__accession__in", accessions,
                             "targets", chunkSize = chunkSize, verbose = verbose)
        rows <- lapply(recs, function(t) {
            comps <- t$target_components %||% list()
            compAcc <- vapply(comps, function(c) c$accession %||% NA_character_,
                              character(1))
            matched <- intersect(compAcc, accessions)
            if (length(matched) == 0L) return(NULL)
            do.call(rbind, lapply(matched, function(acc) {
                comp <- Filter(function(c) identical(c$accession %||% NA, acc),
                               comps)
                comp <- if (length(comp)) comp[[1]] else list()
                data.frame(
                    QueryIDs         = acc,
                    ChEMBL_TID       = t$target_chembl_id %||% NA_character_,
                    Organism         = t$organism %||% NA_character_,
                    UniProt_ID       = acc,
                    Desc             = comp$component_description %||%
                                       NA_character_,
                    stringsAsFactors = FALSE)
            }))
        })
        out <- do.call(rbind, rows)
        if (is.null(out)) {
            out <- data.frame(QueryIDs = character(0), ChEMBL_TID = character(0),
                Organism = character(0), UniProt_ID = character(0),
                Desc = character(0), stringsAsFactors = FALSE)
        }
        out
    }

    ## Resolve ChEMBL_TID -> UniProt accession(s)/Organism/Desc (drug ->
    ## target direction, after mechanism lookup has produced target IDs).
    .resolveTargetsById <- function(targetChemblIds) {
        empty <- data.frame(ChEMBL_TID = character(0), UniProt_ID = character(0),
            Organism = character(0), Desc = character(0), stringsAsFactors = FALSE)
        recs <- .dtiBatchGET(paste0(base, "/target.json"), "target_chembl_id__in",
                             targetChemblIds, "targets", chunkSize = chunkSize,
                             verbose = verbose)
        if (length(recs) == 0L) return(empty)
        rows <- lapply(recs, function(t) {
            comps <- t$target_components %||% list()
            if (length(comps) == 0L) {
                return(data.frame(
                    ChEMBL_TID = t$target_chembl_id %||% NA_character_,
                    UniProt_ID = NA_character_,
                    Organism = t$organism %||% NA_character_,
                    Desc = NA_character_, stringsAsFactors = FALSE))
            }
            do.call(rbind, lapply(comps, function(c) data.frame(
                ChEMBL_TID = t$target_chembl_id %||% NA_character_,
                UniProt_ID = c$accession %||% NA_character_,
                Organism   = t$organism %||% NA_character_,
                Desc       = c$component_description %||% NA_character_,
                stringsAsFactors = FALSE)))
        })
        do.call(rbind, rows)
    }

    ## Mesh indication terms per molecule, collapsed the same way the local
    ## SQL query's GROUP_CONCAT(mesh_id || ': ' || mesh_heading) does.
    .meshFor <- function(moleculeChemblIds) {
        moleculeChemblIds <- unique(stats::na.omit(moleculeChemblIds))
        out <- stats::setNames(rep(NA_character_, length(moleculeChemblIds)),
                        moleculeChemblIds)
        recs <- .dtiBatchGET(paste0(base, "/drug_indication.json"),
                             "molecule_chembl_id__in", moleculeChemblIds,
                             "drug_indications", chunkSize = chunkSize,
                             verbose = verbose)
        if (length(recs) == 0L) return(out)
        df <- do.call(rbind, lapply(recs, function(i) data.frame(
            molecule_chembl_id = i$molecule_chembl_id %||% NA_character_,
            mesh = if (!is.null(i$mesh_id) && !is.null(i$mesh_heading))
                paste0(i$mesh_id, ": ", i$mesh_heading) else NA_character_,
            stringsAsFactors = FALSE)))
        agg <- tapply(df$mesh, df$molecule_chembl_id,
                     function(x) paste(unique(stats::na.omit(x)), collapse = "; "))
        out[names(agg)] <- agg
        out
    }

    ## drug_mechanism equivalent. Batched via <field>__in; the returned
    ## records already carry target_chembl_id/molecule_chembl_id, which
    ## doubles as the join/query key, so no separate per-ID query tag is
    ## needed the way a one-request-per-ID loop would require.
    .mechanismsFor <- function(param, ids) {
        recs <- .dtiBatchGET(paste0(base, "/mechanism.json"), paste0(param, "__in"),
                             ids, "mechanisms", chunkSize = chunkSize,
                             verbose = verbose)
        if (length(recs) == 0L) return(NULL)
        do.call(rbind, lapply(recs, function(m) data.frame(
            chembl_id        = m$molecule_chembl_id %||% NA_character_,
            ChEMBL_TID       = m$target_chembl_id   %||% NA_character_,
            MOA              = m$mechanism_of_action %||% NA_character_,
            Action_Type      = m$action_type        %||% NA_character_,
            Max_Phase        = as.numeric(m$max_phase %||% NA),
            stringsAsFactors = FALSE)))
    }

    if (isTarget) {
        tgtByAcc <- .resolveTargetsByAccession(queryBy$ids)
        if (nrow(tgtByAcc) == 0L) return(emptyDF)
        mech <- .mechanismsFor("target_chembl_id", unique(tgtByAcc$ChEMBL_TID))
        if (is.null(mech) || nrow(mech) == 0L) return(emptyDF)
        out <- merge(tgtByAcc, mech, by = "ChEMBL_TID")
    } else {
        mech <- .mechanismsFor("molecule_chembl_id", queryBy$ids)
        if (is.null(mech) || nrow(mech) == 0L) return(emptyDF)
        mech$QueryIDs <- mech$chembl_id
        tgtById <- .resolveTargetsById(unique(mech$ChEMBL_TID))
        out <- merge(mech, tgtById, by = "ChEMBL_TID", all.x = TRUE)
    }

    drugNames <- getChemblMolecule(unique(out$chembl_id), chunkSize = chunkSize)
    out <- merge(out, drugNames[, c("chembl_id", "pref_name", "first_approval")],
                by = "chembl_id", all.x = TRUE)
    mesh <- .meshFor(unique(out$chembl_id))
    out$Mesh_Indication <- unname(mesh[out$chembl_id])
    names(out)[names(out) == "pref_name"] <- "Drug_Name"
    names(out)[names(out) == "first_approval"] <- "First_Approval"

    ## Unmatched query IDs still surface as an all-NA row, rather than
    ## silently vanishing, mirroring drugTargetAnnot()'s QueryIDs handling.
    unmatched <- setdiff(queryBy$ids, unique(out$QueryIDs))
    if (length(unmatched)) {
        extra <- out[rep(NA_integer_, length(unmatched)), , drop = FALSE]
        extra$QueryIDs <- unmatched
        out <- rbind(out, extra)
    }
    out <- out[order(match(out$QueryIDs, queryBy$ids)), emptyCols]
    rownames(out) <- NULL
    out
}


## ---------------------------------------------------------------------
## PubChem PUG-REST + NCBI E-utilities access
## ---------------------------------------------------------------------
## PubChem's bioactivity data lives in one large "concise" table per
## query entity: gene-centric (/gene/geneid/<id>/concise/JSON) for
## target -> drug, or compound-centric (/compound/cid/<cid>/assaysummary/
## JSON) for drug -> target. Both share the same column shape (Activity
## Outcome, Activity Name, "Activity Value [uM]", Target Accession,
## Target GeneID, CID, Assay Name) and the same Active/numeric/potency
## filtering logic, ported from the validated R_Py_code/pubchem_fetch.py
## reference (gene-centric direction only).
##
## Two things PubChem itself doesn't resolve, confirmed live 2026-07-13:
##  1. Gene-centric queries need a human NCBI GeneID up front - the
##     PUG-REST /gene/symbol route is unreliable for this, so symbols are
##     resolved via E-utilities esearch restricted to one taxon (mirrors
##     the Python reference exactly).
##  2. Compound-centric queries return GeneIDs for WHATEVER species each
##     assay used - e.g. aspirin's classic COX1/COX2 potency assays are
##     annotated against *Ovis aries* (sheep) GeneIDs, not human
##     PTGS1/PTGS2. Those GeneIDs are resolved back to (symbol, taxid)
##     via E-utilities esummary and filtered to `taxid` (default 9606 =
##     human), matching the gene-centric side's scope.
##
## Unlike ChEMBL, PubChem's target and compound directions are genuinely
## asymmetric (no shared `<field>__in=` filter covers both), so - as with
## Open Targets - target->drug and drug->target stay separate functions
## (getPubchemDrugs / getPubchemTargets) rather than one dispatcher, per
## the drug_mechanism-style join ChEMBL supports. getPubchemDrugTarget()
## below adds a thin queryBy-uniform wrapper over both for interface
## consistency with getChemblDrugTarget() / drugTargetAnnot().

.dtiPotencyEndpoints <- c("IC50", "Ki", "Kd", "EC50", "AC50", "Potency")

#' Resolve a gene symbol to an NCBI GeneID via E-utilities esearch
#'
#' The PUG-REST \code{/gene/symbol/...} route is unreliable for this
#' (per the validated Python reference); esearch restricted to one taxon
#' is the robust path.
#'
#' @param symbol character(1) HGNC gene symbol.
#' @param taxid integer(1) NCBI taxonomy ID (default 9606 = human).
#' @return character(1) NCBI GeneID, or \code{NA} if unresolved.
#' @keywords internal
.dtiPubchemGeneId <- function(symbol, taxid = 9606L) {
    url <- paste0(.dtiEndpoints()$eutils, "/esearch.fcgi")
    res <- .dtiApiGET(url, query = list(
        db = "gene", term = sprintf("%s[sym] AND %d[taxid]", symbol, taxid),
        retmode = "json"))
    ids <- res$esearchresult$idlist %||% list()
    if (length(ids) == 0L) return(NA_character_)
    as.character(ids[[1]])
}

#' Raw concise bioactivity table for one GeneID
#' @keywords internal
.dtiPubchemGeneConcise <- function(geneid) {
    url <- paste0(.dtiEndpoints()$pubchem, "/gene/geneid/", geneid, "/concise/JSON")
    res <- .dtiApiGET(url)
    tbl <- res$Table %||% list()
    list(cols = tbl$Columns$Column %||% list(), rows = tbl$Row %||% list())
}

#' Raw compound-centric assay-summary table for one or more CIDs (batched)
#' @keywords internal
.dtiPubchemAssaySummary <- function(cids, chunkSize = 50L) {
    base <- .dtiEndpoints()$pubchem
    cols <- NULL; allRows <- list()
    for (ch in .dtiChunk(cids, chunkSize)) {
        url <- paste0(base, "/compound/cid/", paste(ch, collapse = ","), "/assaysummary/JSON")
        res <- .dtiApiGET(url)
        tbl <- res$Table %||% list()
        if (is.null(cols)) cols <- tbl$Columns$Column %||% list()
        allRows <- c(allRows, tbl$Row %||% list())
    }
    list(cols = cols, rows = allRows)
}

#' Filter+tidy raw concise/assaysummary rows to Active, numeric, potency
#'
#' Shared by both directions - PubChem's gene-centric "concise" table and
#' compound-centric "assaysummary" table use the same column names.
#'
#' @param cols list of raw column-name strings (PUG-REST \code{Columns}).
#' @param rows list of raw row objects (PUG-REST \code{Row}, each with a
#'   \code{Cell} list).
#' @param keepActivities character vector of Activity Name values to
#'   keep, or \code{NULL} to keep all.
#' @return data.frame with columns \code{CID}, \code{TargetAccession},
#'   \code{TargetGeneID}, \code{ActivityName}, \code{ActivityValueuM},
#'   \code{AssayName}.
#' @keywords internal
.dtiPubchemFilterRows <- function(cols, rows, keepActivities = .dtiPotencyEndpoints) {
    cols <- vapply(cols, function(x) x, character(1))
    idx <- stats::setNames(seq_along(cols), cols)
    get <- function(cell, name) {
        i <- idx[name]
        if (is.na(i) || i > length(cell)) return(NA_character_)
        v <- cell[[i]]
        if (is.null(v) || identical(v, "")) NA_character_ else as.character(v)
    }
    recs <- lapply(rows, function(rw) {
        cell <- rw$Cell %||% list()
        outcome <- get(cell, "Activity Outcome")
        if (is.na(outcome) || tolower(outcome) != "active") return(NULL)
        aname <- get(cell, "Activity Name")
        if (!is.null(keepActivities) && (is.na(aname) || !(aname %in% keepActivities))) return(NULL)
        val <- suppressWarnings(as.numeric(get(cell, "Activity Value [uM]")))
        if (is.na(val)) return(NULL)
        data.frame(
            CID              = get(cell, "CID"),
            TargetAccession  = get(cell, "Target Accession"),
            TargetGeneID     = get(cell, "Target GeneID"),
            ActivityName     = aname,
            ActivityValueuM  = val,
            AssayName        = get(cell, "Assay Name"),
            stringsAsFactors = FALSE)
    })
    recs <- recs[!vapply(recs, is.null, logical(1))]
    if (length(recs) == 0L) {
        return(data.frame(CID = character(), TargetAccession = character(),
                          TargetGeneID = character(), ActivityName = character(),
                          ActivityValueuM = numeric(), AssayName = character(),
                          stringsAsFactors = FALSE))
    }
    do.call(rbind, recs)
}

#' Keep the most potent (lowest activity value) row per CID, capped
#' @keywords internal
.dtiPubchemMostPotentPerCid <- function(df, maxCids = NULL) {
    df <- df[order(df$ActivityValueuM), , drop = FALSE]
    df <- df[!duplicated(df$CID), , drop = FALSE]
    if (!is.null(maxCids) && nrow(df) > maxCids) df <- df[seq_len(maxCids), , drop = FALSE]
    df
}

#' CID -> compound Title (drug/common name) + SMILES, batched
#' @keywords internal
.dtiPubchemCidProps <- function(cids, chunkSize = 100L) {
    base <- .dtiEndpoints()$pubchem
    cids <- unique(cids[!is.na(cids) & nzchar(cids)])
    empty <- data.frame(CID = character(), Title = character(),
                        SMILES = character(), stringsAsFactors = FALSE)
    if (length(cids) == 0L) return(empty)
    accum <- list()
    for (ch in .dtiChunk(cids, chunkSize)) {
        url <- paste0(base, "/compound/cid/", paste(ch, collapse = ","),
                      "/property/Title,SMILES/JSON")
        res <- .dtiApiGET(url)
        props <- res$PropertyTable$Properties %||% list()
        if (length(props) > 0L) accum <- c(accum, lapply(props, function(p) data.frame(
            CID    = as.character(p$CID %||% NA_character_),
            Title  = p$Title  %||% NA_character_,
            SMILES = p$SMILES %||% NA_character_,
            stringsAsFactors = FALSE)))
    }
    if (length(accum) == 0L) return(empty)
    do.call(rbind, accum)
}

#' Resolve a drug/compound name to a PubChem CID
#' @param name character(1) compound name, e.g. \code{"aspirin"}.
#' @return character(1) CID, or \code{NA} if unresolved.
#' @keywords internal
.dtiPubchemNameToCid <- function(name) {
    url <- paste0(.dtiEndpoints()$pubchem, "/compound/name/",
                  utils::URLencode(name, reserved = TRUE), "/cids/JSON")
    res <- .dtiApiGET(url)
    ids <- res$IdentifierList$CID %||% list()
    if (length(ids) == 0L) return(NA_character_)
    as.character(ids[[1]])
}

#' GeneID -> (symbol, taxid) lookup via E-utilities esummary, batched
#' @keywords internal
.dtiPubchemGeneInfo <- function(geneids) {
    empty <- data.frame(GeneID = character(), Symbol = character(),
                        Taxid = integer(), stringsAsFactors = FALSE)
    geneids <- unique(geneids[!is.na(geneids) & nzchar(geneids)])
    if (length(geneids) == 0L) return(empty)
    url <- paste0(.dtiEndpoints()$eutils, "/esummary.fcgi")
    accum <- list()
    for (ch in .dtiChunk(geneids, 200L)) {
        res <- .dtiApiGET(url, query = list(db = "gene",
                                            id = paste(ch, collapse = ","),
                                            retmode = "json"))
        result <- res$result %||% list()
        uids <- result$uids %||% list()
        for (u in uids) {
            g <- result[[as.character(u)]]
            accum[[length(accum) + 1L]] <- data.frame(
                GeneID = as.character(u),
                Symbol = g$name %||% NA_character_,
                Taxid  = as.integer(g$organism$taxid %||% NA),
                stringsAsFactors = FALSE)
        }
    }
    if (length(accum) == 0L) return(empty)
    do.call(rbind, accum)
}

#' Retrieve PubChem bioactivities for one or more genes (target -> drug)
#'
#' Ports the validated \code{R_Py_code/pubchem_fetch.py} reference:
#' resolves each gene symbol to a human NCBI GeneID, pulls its concise
#' bioactivity table, keeps Active + numeric + potency-endpoint rows,
#' collapses to the most potent row per CID (capped at \code{maxCids}),
#' and looks up each surviving CID's name and SMILES.
#'
#' @param genes character vector of HGNC gene symbols (or NCBI GeneIDs,
#'   which pass through unresolved).
#' @param taxid integer(1) NCBI taxonomy ID for symbol resolution
#'   (default 9606 = human).
#' @param keepActivities character vector of Activity Name values to
#'   keep (default IC50/Ki/Kd/EC50/AC50/Potency); \code{NULL} keeps all.
#' @param maxCids integer(1) cap on CIDs kept per gene (default 400).
#' @param pause numeric(1) seconds between requests (NCBI: <=3 req/s
#'   without an API key).
#' @param verbose logical(1) progress messages.
#' @return A \code{data.frame} with columns \code{gene_symbol},
#'   \code{geneid}, \code{target_accession}, \code{cid}, \code{drug_name}
#'   (~27% are patent references, not common names - inherent to PubChem
#'   breadth), \code{canonical_smiles}, \code{activity_name},
#'   \code{activity_value_uM}, \code{assay_name}, \code{db}. Empty if
#'   none / offline.
#' @examples
#' \donttest{
#'   df <- getPubchemDrugs(c("FGFR1", "KLB"))
#'   table(df$gene_symbol)
#' }
#' @seealso \code{\link{getPubchemTargets}}, \code{\link{getPubchemDrugTarget}}
#' @export
getPubchemDrugs <- function(genes, taxid = 9606L,
                            keepActivities = .dtiPotencyEndpoints,
                            maxCids = 400L, pause = 0.34, verbose = FALSE) {
    stopifnot(is.character(genes), length(genes) >= 1L)
    rows <- vector("list", length(genes))
    for (i in seq_along(genes)) {
        g <- genes[i]
        gid <- if (grepl("^[0-9]+$", g)) g else .dtiPubchemGeneId(g, taxid = taxid)
        if (verbose) message("PubChem gene ", g, " -> GeneID ", gid)
        if (!is.na(gid)) {
            raw <- .dtiPubchemGeneConcise(gid)
            df <- .dtiPubchemFilterRows(raw$cols, raw$rows, keepActivities)
            if (nrow(df) > 0L) {
                df <- .dtiPubchemMostPotentPerCid(df, maxCids)
                props <- .dtiPubchemCidProps(df$CID)
                df$drug_name         <- props$Title[match(df$CID, props$CID)]
                df$canonical_smiles  <- props$SMILES[match(df$CID, props$CID)]
                df$gene_symbol       <- g
                df$geneid            <- gid
                rows[[i]] <- df
            }
        }
        Sys.sleep(pause)
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) {
        return(as.data.frame(stats::setNames(replicate(10, character(0), simplify = FALSE),
            c("gene_symbol", "geneid", "target_accession", "cid", "drug_name",
              "canonical_smiles", "activity_name", "activity_value_uM",
              "assay_name", "db")), stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, rows)
    out <- out[, c("gene_symbol", "geneid", "TargetAccession", "CID", "drug_name",
                   "canonical_smiles", "ActivityName", "ActivityValueuM", "AssayName")]
    colnames(out) <- c("gene_symbol", "geneid", "target_accession", "cid", "drug_name",
                       "canonical_smiles", "activity_name", "activity_value_uM", "assay_name")
    out$db <- "PubChem"
    rownames(out) <- NULL
    out
}

#' Retrieve PubChem bioactivities for one or more drugs (drug -> target)
#'
#' Drug -> target counterpart of \code{\link{getPubchemDrugs}}: resolves
#' each drug name to a CID, pulls its compound-centric assay-summary
#' table (batched), keeps Active + numeric + potency-endpoint rows,
#' collapses to the most potent row per (CID, GeneID) pair, and resolves
#' each target GeneID to a gene symbol - filtered to \code{taxid}
#' (default human), since compound assay data spans whatever species
#' each assay used (see file header note).
#'
#' @param drugs character vector of drug/compound names (or PubChem
#'   CIDs, which pass through unresolved).
#' @param taxid integer(1) or \code{NULL}; keep only targets from this
#'   NCBI taxonomy ID (default 9606 = human), or all species if
#'   \code{NULL}.
#' @param keepActivities character vector of Activity Name values to
#'   keep (default IC50/Ki/Kd/EC50/AC50/Potency); \code{NULL} keeps all.
#' @param maxCidsPerBatch integer(1) CIDs per \code{assaysummary} HTTP
#'   request (default 50).
#' @param pause numeric(1) seconds between name-resolution requests.
#' @param verbose logical(1) progress messages.
#' @return A \code{data.frame} with columns \code{cid}, \code{drug_name},
#'   \code{gene_symbol}, \code{geneid}, \code{taxid},
#'   \code{target_accession}, \code{activity_name},
#'   \code{activity_value_uM}, \code{assay_name}, \code{db}. Empty if
#'   none / offline.
#' @examples
#' \donttest{
#'   df <- getPubchemTargets("aspirin")
#'   df[, c("drug_name", "gene_symbol", "activity_name", "activity_value_uM")]
#' }
#' @seealso \code{\link{getPubchemDrugs}}, \code{\link{getPubchemDrugTarget}}
#' @export
getPubchemTargets <- function(drugs, taxid = 9606L,
                              keepActivities = .dtiPotencyEndpoints,
                              maxCidsPerBatch = 50L, pause = 0.2, verbose = FALSE) {
    stopifnot(is.character(drugs), length(drugs) >= 1L)
    empty <- as.data.frame(stats::setNames(replicate(10, character(0), simplify = FALSE),
        c("cid", "drug_name", "gene_symbol", "geneid", "taxid",
          "target_accession", "activity_name", "activity_value_uM",
          "assay_name", "db")), stringsAsFactors = FALSE)

    cidFor <- character(length(drugs))
    for (i in seq_along(drugs)) {
        d <- drugs[i]
        cidFor[i] <- if (grepl("^[0-9]+$", d)) d else .dtiPubchemNameToCid(d)
        if (verbose) message("PubChem drug ", d, " -> CID ", cidFor[i])
        Sys.sleep(pause)
    }
    names(cidFor) <- drugs
    cids <- unique(stats::na.omit(cidFor))
    if (length(cids) == 0L) return(empty)

    raw <- .dtiPubchemAssaySummary(cids, chunkSize = maxCidsPerBatch)
    df  <- .dtiPubchemFilterRows(raw$cols, raw$rows, keepActivities)
    if (nrow(df) == 0L) return(empty)
    df <- df[order(df$ActivityValueuM), , drop = FALSE]
    df <- df[!duplicated(paste(df$CID, df$TargetGeneID)), , drop = FALSE]

    genes <- .dtiPubchemGeneInfo(df$TargetGeneID)
    df$gene_symbol <- genes$Symbol[match(df$TargetGeneID, genes$GeneID)]
    df$taxid       <- genes$Taxid[match(df$TargetGeneID, genes$GeneID)]
    if (!is.null(taxid)) df <- df[!is.na(df$taxid) & df$taxid == taxid, , drop = FALSE]
    if (nrow(df) == 0L) return(empty)

    drugNameFor <- stats::setNames(names(cidFor), cidFor)
    df$drug_name <- unname(drugNameFor[df$CID])

    out <- df[, c("CID", "drug_name", "gene_symbol", "TargetGeneID", "taxid",
                 "TargetAccession", "ActivityName", "ActivityValueuM", "AssayName")]
    colnames(out) <- c("cid", "drug_name", "gene_symbol", "geneid", "taxid",
                       "target_accession", "activity_name", "activity_value_uM", "assay_name")
    out$db <- "PubChem"
    rownames(out) <- NULL
    out
}

#' Query PubChem bioactivity data via the uniform queryBy interface
#'
#' Thin \code{queryBy}-dispatching wrapper over \code{\link{getPubchemDrugs}}
#' (target -> drug) and \code{\link{getPubchemTargets}} (drug -> target),
#' matching the \code{queryBy = list(molType, idType, ids)} interface used by
#' \code{\link{drugTargetAnnot}} and \code{\link{getChemblDrugTarget}}, so a
#' future cross-source ID-translation/meta layer can dispatch to any source
#' function the same way. PubChem's two directions are genuinely asymmetric
#' (no shared bulk filter covers both, unlike ChEMBL's
#' \code{drug_mechanism} join), so this only adds a \code{QueryIDs} column
#' and unmatched-ID NA-row padding on top of the two underlying functions -
#' it does not change their query logic.
#'
#' @param queryBy named list with components \code{molType}, \code{idType}
#'   and \code{ids}. \code{molType = "gene"} with \code{idType = "symbol"}
#'   queries target -> drug; \code{ids} may be HGNC gene symbols or NCBI
#'   GeneIDs (numeric strings pass through unresolved), matching
#'   \code{\link{getPubchemDrugs}}. \code{molType = "cmp"} with
#'   \code{idType = "name"} queries drug -> target; \code{ids} may be
#'   compound names or PubChem CIDs (numeric strings pass through
#'   unresolved), matching \code{\link{getPubchemTargets}}. Other identifier
#'   types (UniProt accession, ChEMBL ID, PubChem CID as the *only* form,
#'   DrugBank ID, ...) are not yet supported over this wrapper.
#' @param ... additional arguments passed through to
#'   \code{\link{getPubchemDrugs}} / \code{\link{getPubchemTargets}} (e.g.
#'   \code{taxid}, \code{keepActivities}, \code{verbose}).
#' @return A \code{data.frame} in the same column shape as
#'   \code{\link{getPubchemDrugs}} / \code{\link{getPubchemTargets}}, plus a
#'   leading \code{QueryIDs} column echoing the original \code{queryBy$ids}
#'   token each row resolved from. Query IDs that returned no rows still
#'   appear as a single row with all other fields \code{NA}, so callers can
#'   always confirm which of their input IDs were resolved.
#' @examples
#' \donttest{
#'   ## target -> drug: FGFR1, KLB
#'   getPubchemDrugTarget(list(molType = "gene", idType = "symbol",
#'                             ids = c("FGFR1", "KLB")))
#'   ## drug -> target: aspirin
#'   getPubchemDrugTarget(list(molType = "cmp", idType = "name",
#'                             ids = "aspirin"))
#' }
#' @seealso \code{\link{getPubchemDrugs}}, \code{\link{getPubchemTargets}},
#'   \code{\link{getChemblDrugTarget}}
#' @export
getPubchemDrugTarget <- function(queryBy = list(molType = NULL, idType = NULL,
                                                ids = NULL), ...) {
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

    isGene <- identical(queryBy$molType, "gene") &&
        identical(queryBy$idType, "symbol")
    isCmp <- identical(queryBy$molType, "cmp") &&
        identical(queryBy$idType, "name")
    if (!isGene && !isCmp) {
        stop(
            "getPubchemDrugTarget() currently supports only ",
            "queryBy=list(molType=\"gene\", idType=\"symbol\", ids=...) ",
            "or queryBy=list(molType=\"cmp\", idType=\"name\", ids=...). ",
            "Other identifier types require translating to a gene symbol ",
            "or a compound name/CID first."
        )
    }

    ids <- queryBy$ids
    if (isGene) {
        out <- getPubchemDrugs(ids, ...)
        queryCol <- "gene_symbol"
    } else {
        out <- getPubchemTargets(ids, ...)
        queryCol <- "drug_name"
    }
    out$QueryIDs <- out[[queryCol]]

    unmatched <- setdiff(ids, unique(out$QueryIDs))
    if (length(unmatched)) {
        extra <- out[rep(NA_integer_, length(unmatched)), , drop = FALSE]
        extra$QueryIDs <- unmatched
        out <- rbind(out, extra)
    }
    front <- c("QueryIDs", setdiff(names(out), "QueryIDs"))
    out <- out[order(match(out$QueryIDs, ids)), front]
    rownames(out) <- NULL
    out
}
