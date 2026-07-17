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

#' Perform a GraphQL POST with retry + polite throttling
#'
#' Shared by any GraphQL-backed source (DGIdb, and later Open Targets).
#' Mirrors \code{\link{.dtiApiGET}}'s retry/throttle/timeout behavior but
#' posts a \code{\{query, variables\}} JSON body instead of a query-string
#' GET. On any transport failure or a GraphQL \code{errors} payload it
#' warns and returns \code{NULL} (or the partial \code{data} field, for
#' errors) so callers can degrade gracefully, unless \code{hardStop = TRUE}.
#'
#' @param url character(1) GraphQL endpoint URL.
#' @param query character(1) GraphQL query/mutation document.
#' @param variables named list of GraphQL variables.
#' @param timeout numeric(1) per-request timeout in seconds.
#' @param maxTries integer(1) total attempts including the first.
#' @param hardStop logical(1) if TRUE, rethrow the error instead of NULL.
#' @return parsed \code{data} element of the GraphQL response, or NULL.
#' @keywords internal
.dtiGraphQL <- function(url, query, variables = list(), timeout = 60L,
                        maxTries = 3L, hardStop = FALSE) {
    out <- tryCatch({
        body <- list(query = query, variables = variables)
        req <- httr2::request(url)
        req <- httr2::req_headers(req, Accept = "application/json",
                                  `Content-Type` = "application/json")
        req <- httr2::req_user_agent(
            req, "drugTargetInteractions R package (Bioconductor)")
        req <- httr2::req_timeout(req, timeout)
        req <- httr2::req_body_json(req, body)
        req <- httr2::req_throttle(req, rate = 5, fill_time_s = 1)
        req <- httr2::req_retry(
            req, max_tries = maxTries,
            is_transient = function(resp)
                httr2::resp_status(resp) %in% c(429L, 500L, 502L, 503L, 504L))
        resp <- httr2::req_perform(req)
        parsed <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        if (!is.null(parsed$errors)) {
            msgs <- vapply(parsed$errors, function(e) e$message %||% NA_character_,
                          character(1))
            warning("drugTargetInteractions GraphQL query to '", url,
                    "' returned errors: ",
                    paste(stats::na.omit(msgs), collapse = "; "), call. = FALSE)
        }
        parsed$data
    }, error = function(e) {
        if (hardStop) stop(e)
        warning("drugTargetInteractions GraphQL POST failed for '", url, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
    })
    out
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


## ---------------------------------------------------------------------
## DGIdb GraphQL access
## ---------------------------------------------------------------------
## Unlike ChEMBL/PubChem, DGIdb's schema is bidirectional by construction:
## a single root query, interactions(geneNames: [...], drugNames: [...]),
## returns Interaction nodes that each already carry BOTH endpoints'
## identity (gene AND drug) plus interactionScore, evidenceScore,
## interactionTypes and sources. So one shared query+parser covers both
## target -> drug and drug -> target - no separate ID-resolution step is
## needed the way ChEMBL/PubChem require (geneNames/drugNames match
## case-insensitively on plain names directly; confirmed live against the
## v5 API 2026-07-13, incl. an exact 227-row match for the project's 8
## genes against the validated R_Py_code/dgidb_fetch.py reference).
##
## DGIdb normalizes matched names to its own canonical casing regardless
## of query casing (confirmed live 2026-07-16: querying "fgfr1" returns
## gene_name "FGFR1"; "Aspirin" returns drug_name "ASPIRIN") - so mapping
## a result row back to the query token it came from (for
## getDgidbDrugTarget()'s QueryIDs column) has to match
## case-insensitively, unlike ChEMBL/PubChem where the query token is
## always echoed back verbatim.
##
## The `sources` column (upstream sourceDbName per interaction, e.g.
## "ChEMBL", "DrugBank", "TTD") is retained specifically because DGIdb is
## an aggregator - individual rows carry their upstream source's license,
## not one blanket DGIdb license (see the standalone
## R_Py_code/license_registry.py table; not ported into the package).

#' Column order for tidy DGIdb interaction rows
#' @keywords internal
.dtiDgidbCols <- c("gene_name", "drug_name", "drug_concept_id", "drug_approved",
                   "interaction_types", "directionality", "interaction_score",
                   "evidence_score", "sources", "db")

#' Empty DGIdb data.frame with the canonical columns
#' @keywords internal
.dtiEmptyDgidb <- function() {
    as.data.frame(stats::setNames(
        replicate(length(.dtiDgidbCols), character(0), simplify = FALSE),
        .dtiDgidbCols), stringsAsFactors = FALSE)
}

#' The reusable GraphQL selection for one Interaction node (no outer braces)
#' @keywords internal
.dtiDgidbSelection <- function() {
    "drug { name conceptId approved }
     gene { name }
     interactionScore
     evidenceScore
     interactionTypes { type directionality }
     sources { sourceDbName }"
}

#' Parse a list of Interaction nodes into a tidy data.frame
#' @keywords internal
.dtiParseDgidbInteractions <- function(nodes) {
    nodes <- nodes %||% list()
    if (length(nodes) == 0L) return(.dtiEmptyDgidb())
    collapse <- function(x) paste(unique(stats::na.omit(x)), collapse = "; ")
    rows <- lapply(nodes, function(n) {
        drug   <- n$drug %||% list()
        gene   <- n$gene %||% list()
        itypes <- n$interactionTypes %||% list()
        types  <- collapse(vapply(itypes, function(t) t$type %||% NA_character_,
                                  character(1)))
        dirs   <- collapse(vapply(itypes, function(t) t$directionality %||% NA_character_,
                                  character(1)))
        srcs   <- n$sources %||% list()
        srcNames <- collapse(vapply(srcs, function(s) s$sourceDbName %||% NA_character_,
                                    character(1)))
        data.frame(
            gene_name         = gene$name %||% NA_character_,
            drug_name         = drug$name %||% NA_character_,
            drug_concept_id   = drug$conceptId %||% NA_character_,
            drug_approved     = as.logical(drug$approved %||% NA),
            interaction_types = types,
            directionality    = dirs,
            interaction_score = as.numeric(n$interactionScore %||% NA),
            evidence_score    = as.numeric(n$evidenceScore %||% NA),
            sources           = srcNames,
            db                = "DGIdb",
            stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    df
}

#' Paginated, chunked interactions() query, shared by both direction wrappers
#'
#' \code{names} is sent as a single GraphQL array variable per request (no
#' URL-length constraint the way ChEMBL's REST \code{__in} filters have),
#' but is still chunked defensively at \code{chunkSize} - consistent with
#' the batching convention used for the other REST sources - since DGIdb
#' publishes no documented cap on array size. Each chunk is paginated via
#' its own cursor until exhausted.
#'
#' @param names character vector of gene symbols or drug names.
#' @param by character(1) \code{"gene"} or \code{"drug"} - which
#'   \code{interactions()} filter argument \code{names} populates.
#' @param pageSize integer(1) rows per GraphQL page (cursor-paginated).
#' @param maxRows integer(1) cap on total returned rows.
#' @param chunkSize integer(1) max names sent per request.
#' @param verbose logical(1); if TRUE, message progress per chunk/page.
#' @return list of raw \code{Interaction} GraphQL nodes (unparsed).
#' @keywords internal
.dtiDgidbFetch <- function(names, by = c("gene", "drug"), pageSize = 500L,
                           maxRows = 5000L, chunkSize = 300L, verbose = FALSE) {
    by <- match.arg(by)
    url <- .dtiEndpoints()$dgidb
    argName <- if (by == "gene") "geneNames" else "drugNames"
    query <- sprintf("
      query($names: [String!], $after: String, $first: Int) {
        interactions(%s: $names, after: $after, first: $first) {
          pageInfo { hasNextPage endCursor }
          nodes { %s }
        }
      }", argName, .dtiDgidbSelection())

    acc <- list()
    for (chunk in .dtiChunk(unique(names), chunkSize)) {
        got <- 0L; after <- NULL
        repeat {
            vars <- list(names = as.list(chunk), first = pageSize, after = after)
            if (verbose)
                message("DGIdb ", argName, " chunk n=", length(chunk),
                        " after=", after %||% "<start>")
            data <- .dtiGraphQL(url, query, variables = vars)
            conn <- data$interactions
            nodes <- conn$nodes %||% list()
            if (length(nodes) == 0L) break
            acc <- c(acc, nodes)
            got <- got + length(nodes)
            hasNext <- isTRUE(conn$pageInfo$hasNextPage %||% FALSE)
            after <- conn$pageInfo$endCursor %||% NULL
            if (!hasNext || got >= maxRows || is.null(after)) break
        }
    }
    acc
}

#' Retrieve DGIdb drug interactions for one or more genes
#'
#' Queries the DGIdb \code{interactions(geneNames: ...)} GraphQL field
#' and returns a tidy \code{data.frame} of gene-drug interaction rows.
#' This is the target -> drug direction; see \code{\link{getDgidbTargets}}
#' for the reverse.
#'
#' @param genes character vector of HGNC gene symbols (case-insensitive;
#'   unmatched symbols simply contribute no rows, no error).
#' @param pageSize integer(1) rows per GraphQL page (default 500).
#' @param maxRows integer(1) cap on total returned rows (default 5000).
#' @param verbose logical(1); if TRUE, message progress per chunk/page.
#' @return A \code{data.frame} with columns \code{gene_name},
#'   \code{drug_name}, \code{drug_concept_id} (e.g.
#'   \code{"chembl:CHEMBL1201585"}), \code{drug_approved},
#'   \code{interaction_types}, \code{directionality},
#'   \code{interaction_score}, \code{evidence_score}, \code{sources}
#'   (upstream provenance), \code{db}. Empty if none / offline.
#' @examples
#' \donttest{
#'   df <- getDgidbDrugs(c("FGFR1", "KLB"))
#'   table(df$gene_name)
#' }
#' @seealso \code{\link{getDgidbTargets}}, \code{\link{getDgidbDrugTarget}}
#' @export
getDgidbDrugs <- function(genes, pageSize = 500L, maxRows = 5000L, verbose = FALSE) {
    stopifnot(is.character(genes), length(genes) >= 1L)
    nodes <- .dtiDgidbFetch(genes, by = "gene", pageSize = pageSize,
                            maxRows = maxRows, verbose = verbose)
    .dtiParseDgidbInteractions(nodes)
}

#' Retrieve DGIdb target interactions for one or more drugs
#'
#' Queries the DGIdb \code{interactions(drugNames: ...)} GraphQL field
#' and returns a tidy \code{data.frame} of drug-gene interaction rows.
#' This is the drug -> target direction; see \code{\link{getDgidbDrugs}}
#' for the reverse. Same columns/shape as \code{\link{getDgidbDrugs}}
#' since \code{Interaction} nodes always carry both endpoints.
#'
#' @param drugs character vector of drug names (case-insensitive;
#'   unmatched names simply contribute no rows, no error).
#' @param pageSize integer(1) rows per GraphQL page (default 500).
#' @param maxRows integer(1) cap on total returned rows (default 5000).
#' @param verbose logical(1); if TRUE, message progress per chunk/page.
#' @return A \code{data.frame}; see \code{\link{getDgidbDrugs}} for columns.
#' @examples
#' \donttest{
#'   df <- getDgidbTargets(c("imatinib", "aspirin"))
#'   table(df$drug_name, useNA = "no")
#' }
#' @seealso \code{\link{getDgidbDrugs}}, \code{\link{getDgidbDrugTarget}}
#' @export
getDgidbTargets <- function(drugs, pageSize = 500L, maxRows = 5000L, verbose = FALSE) {
    stopifnot(is.character(drugs), length(drugs) >= 1L)
    nodes <- .dtiDgidbFetch(drugs, by = "drug", pageSize = pageSize,
                            maxRows = maxRows, verbose = verbose)
    .dtiParseDgidbInteractions(nodes)
}

#' Query DGIdb interaction data via the uniform queryBy interface
#'
#' Thin \code{queryBy}-dispatching wrapper over \code{\link{getDgidbDrugs}}
#' (target -> drug) and \code{\link{getDgidbTargets}} (drug -> target),
#' matching the \code{queryBy = list(molType, idType, ids)} interface used
#' by \code{\link{drugTargetAnnot}}, \code{\link{getChemblDrugTarget}} and
#' \code{\link{getPubchemDrugTarget}}. Since DGIdb normalizes matched names
#' to its own canonical casing (see the file header note), the returned
#' \code{QueryIDs} column is matched back to the original \code{queryBy$ids}
#' token \emph{case-insensitively} rather than by verbatim string equality.
#'
#' @param queryBy named list with components \code{molType}, \code{idType}
#'   and \code{ids}. \code{molType = "gene"} with \code{idType = "symbol"}
#'   queries target -> drug (\code{ids} = gene symbols). \code{molType =
#'   "cmp"} with \code{idType = "name"} queries drug -> target (\code{ids}
#'   = drug names). Matches \code{\link{getPubchemDrugTarget}}'s
#'   vocabulary for the same concepts.
#' @param ... additional arguments passed through to
#'   \code{\link{getDgidbDrugs}} / \code{\link{getDgidbTargets}} (e.g.
#'   \code{pageSize}, \code{maxRows}, \code{verbose}).
#' @return A \code{data.frame} in the same column shape as
#'   \code{\link{getDgidbDrugs}} / \code{\link{getDgidbTargets}}, plus a
#'   leading \code{QueryIDs} column echoing the original \code{queryBy$ids}
#'   token each row case-insensitively matched to. Query IDs that returned
#'   no rows still appear as a single row with all other fields \code{NA}.
#' @examples
#' \donttest{
#'   ## target -> drug: FGFR1, KLB
#'   getDgidbDrugTarget(list(molType = "gene", idType = "symbol",
#'                           ids = c("FGFR1", "KLB")))
#'   ## drug -> target: imatinib
#'   getDgidbDrugTarget(list(molType = "cmp", idType = "name",
#'                           ids = "imatinib"))
#' }
#' @seealso \code{\link{getDgidbDrugs}}, \code{\link{getDgidbTargets}},
#'   \code{\link{getChemblDrugTarget}}, \code{\link{getPubchemDrugTarget}}
#' @export
getDgidbDrugTarget <- function(queryBy = list(molType = NULL, idType = NULL,
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
            "getDgidbDrugTarget() currently supports only ",
            "queryBy=list(molType=\"gene\", idType=\"symbol\", ids=...) ",
            "or queryBy=list(molType=\"cmp\", idType=\"name\", ids=...). ",
            "Other identifier types require translating to a gene symbol ",
            "or a drug name first."
        )
    }

    ids <- queryBy$ids
    if (isGene) {
        out <- getDgidbDrugs(ids, ...)
        resolvedCol <- "gene_name"
    } else {
        out <- getDgidbTargets(ids, ...)
        resolvedCol <- "drug_name"
    }

    ## DGIdb returns canonically-cased names regardless of query casing
    ## (see file header note), so map each row back to its query token
    ## case-insensitively rather than by verbatim string equality.
    out$QueryIDs <- ids[match(toupper(out[[resolvedCol]]), toupper(ids))]

    unmatched <- ids[!(toupper(ids) %in% toupper(out$QueryIDs))]
    if (length(unmatched)) {
        extra <- out[rep(NA_integer_, length(unmatched)), , drop = FALSE]
        extra$QueryIDs <- unmatched
        out <- rbind(out, extra)
    }
    front <- c("QueryIDs", setdiff(names(out), "QueryIDs"))
    out <- out[order(match(toupper(out$QueryIDs), toupper(ids))), front]
    rownames(out) <- NULL
    out
}


## ---------------------------------------------------------------------
## Open Targets GraphQL access
## ---------------------------------------------------------------------
## Unlike DGIdb, Open Targets' schema is asymmetric: target -> drug goes
## through Target.drugAndClinicalCandidates, drug -> target through
## Drug.mechanismsOfAction.rows[].targets - two separate GraphQL shapes,
## each needing its own ID-resolution step first (gene symbol -> Ensembl
## ID via a `search` query restricted to entity "target"; drug name ->
## ChEMBL ID via `search` restricted to entity "drug" - Open Targets drug
## IDs *are* ChEMBL IDs). Confirmed live 2026-07-17 against a schema
## introspection of the v4 API: Query.drug(chemblId)/Query.drugs(chemblIds)
## are full siblings of Query.target/Query.targets, and
## Drug.mechanismsOfAction.rows[].targets[] can list more than one target
## per mechanism row (e.g. aspirin's "Cyclooxygenase inhibitor" mechanism
## lists both PTGS1 and PTGS2).
##
## Only the proteome-scale/batched shape from the R_Py_code reference is
## ported here (renamed to drop the "Batch" suffix, matching every other
## source's convention of one vector-accepting accessor per direction
## that chunks/aliases internally); the reference's separate single-item
## functions (getOpenTargetsId(), a singular getOpenTargetsDrugs(),
## getOpenTargetsDrugId(), a singular getOpenTargetsTargets()) were
## parallel, non-DRY duplicate implementations of the exact same queries
## for a length-1 input and are not needed once the batched shape handles
## length 1 fine.
##
## getOpenTargetsDrugTarget()'s QueryIDs join uses the *stable resolved
## ID* (ensembl_id / chembl_id), not the display name the way DGIdb's
## wrapper does - Open Targets already returns those IDs verbatim in
## every output row, and joining on an ID sidesteps DGIdb's
## casing-normalization gotcha entirely. This costs one extra batched
## resolution call inside the wrapper (symbols/names it already resolves
## a second time internally) - a deliberate, modest tradeoff for keeping
## getOpenTargetsDrugs()/getOpenTargetsTargets() unchanged in their
## natural, self-contained shape rather than threading extra bookkeeping
## through them for the wrapper's benefit alone.

#' Resolve many gene symbols to Ensembl IDs in batched GraphQL requests
#'
#' Uses GraphQL field-aliasing to resolve up to \code{chunkSize} symbols
#' per HTTP request, so resolving many genes takes a handful of requests
#' rather than one per gene.
#'
#' @param symbols character vector of gene symbols (Ensembl gene IDs pass
#'   through unresolved - checked via a plain \code{"^ENSG[0-9]+$"} regex).
#' @param chunkSize integer(1) symbols per HTTP request (default 100).
#' @param pause numeric(1) seconds to sleep between requests (politeness).
#' @param verbose logical(1) progress messages.
#' @return A named character vector (names = input symbols) of Ensembl
#'   gene IDs; unresolved symbols are \code{NA}. Order matches the input.
#' @examples
#' \donttest{
#'   getOpenTargetsIds(c("FGFR1", "KLB", "NOT_A_GENE"))
#' }
#' @seealso \code{\link{getOpenTargetsDrugs}}
#' @export
getOpenTargetsIds <- function(symbols, chunkSize = 100L, pause = 0.1,
                              verbose = FALSE) {
    stopifnot(is.character(symbols), length(symbols) >= 1L)
    url <- .dtiEndpoints()$opentargets
    uniq <- unique(symbols)
    chunks <- .dtiChunk(uniq, chunkSize)
    resolved <- character(0)
    for (i in seq_along(chunks)) {
        ch <- chunks[[i]]
        if (verbose) message("getOpenTargetsIds chunk ", i, "/", length(chunks),
                             " (", length(ch), " symbols)")
        ## One aliased search per symbol; $qN variables keep it injection-safe.
        varDefs <- paste(sprintf("$q%d: String!", seq_along(ch)), collapse = ", ")
        aliases <- paste(sprintf(
            "a%d: search(queryString: $q%d, entityNames: [\"target\"], page: {index: 0, size: 3}) { hits { id object { ... on Target { approvedSymbol } } } }",
            seq_along(ch), seq_along(ch)), collapse = "\n")
        query <- sprintf("query resolveIds(%s) {\n%s\n}", varDefs, aliases)
        vars <- stats::setNames(as.list(ch), sprintf("q%d", seq_along(ch)))
        data <- .dtiGraphQL(url, query, variables = vars)
        for (j in seq_along(ch)) {
            hits <- data[[sprintf("a%d", j)]]$hits %||% list()
            id <- NA_character_
            if (length(hits) > 0L) {
                exact <- Filter(function(h) {
                    s <- h$object$approvedSymbol %||% NA_character_
                    !is.na(s) && toupper(s) == toupper(ch[j])
                }, hits)
                pick <- if (length(exact) > 0L) exact[[1]] else hits[[1]]
                id <- pick$id %||% NA_character_
            }
            resolved[ch[j]] <- id
        }
        if (pause > 0 && i < length(chunks)) Sys.sleep(pause)
    }
    resolved[symbols]  # re-expand to input order/length
}

#' Column order shared by the drug accessor
#' @keywords internal
.dtiDrugCols <- c("ensembl_id", "approved_symbol", "drug_id", "drug_name",
                  "drug_type", "max_clinical_stage", "mechanism_of_action",
                  "action_type", "disease_id", "disease_name")

#' Empty drug data.frame with the canonical columns
#' @keywords internal
.dtiEmptyDrugs <- function() {
    as.data.frame(stats::setNames(
        replicate(length(.dtiDrugCols), character(0), simplify = FALSE),
        .dtiDrugCols), stringsAsFactors = FALSE)
}

#' The reusable GraphQL selection for one target's drugs (no outer braces)
#' @keywords internal
.dtiDrugSelection <- function() {
    "approvedSymbol
     drugAndClinicalCandidates {
       count
       rows {
         maxClinicalStage
         drug {
           id name drugType maximumClinicalStage
           mechanismsOfAction { rows { mechanismOfAction actionType } }
         }
         diseases { disease { id name } }
       }
     }"
}

#' Parse one target object into a tidy data.frame, applying the requested
#' row expansion.
#'
#' @param tgt list; the GraphQL \code{target} object (may be NULL).
#' @param ensg character(1); the Ensembl ID this object was queried with.
#' @param expand character(1); one of "mechanism", "disease", "drug".
#' @return data.frame with columns \code{.dtiDrugCols}; empty if no rows.
#' @keywords internal
.dtiParseTargetDrugs <- function(tgt, ensg, expand) {
    tgt  <- tgt %||% list()
    rows <- tgt$drugAndClinicalCandidates$rows %||% list()
    if (length(rows) == 0L) return(.dtiEmptyDrugs())
    sym  <- tgt$approvedSymbol %||% NA_character_
    collapse <- function(x) paste(unique(stats::na.omit(x)), collapse = "; ")
    perRow <- lapply(rows, function(r) {
        drug <- r$drug %||% list()
        moaRows <- drug$mechanismsOfAction$rows %||% list()
        moa <- if (length(moaRows) == 0L)
            data.frame(mechanism_of_action = NA_character_,
                       action_type = NA_character_, stringsAsFactors = FALSE)
        else do.call(rbind, lapply(moaRows, function(m) data.frame(
            mechanism_of_action = m$mechanismOfAction %||% NA_character_,
            action_type         = m$actionType        %||% NA_character_,
            stringsAsFactors = FALSE)))
        moa <- unique(moa)
        dis <- r$diseases %||% list()
        dz <- if (length(dis) == 0L)
            data.frame(disease_id = NA_character_,
                       disease_name = NA_character_, stringsAsFactors = FALSE)
        else do.call(rbind, lapply(dis, function(d) data.frame(
            disease_id   = d$disease$id   %||% NA_character_,
            disease_name = d$disease$name %||% NA_character_,
            stringsAsFactors = FALSE)))
        dz <- unique(dz)
        base <- data.frame(
            ensembl_id      = ensg,
            approved_symbol = sym,
            drug_id         = drug$id %||% NA_character_,
            drug_name       = drug$name %||% NA_character_,
            drug_type       = drug$drugType %||% NA_character_,
            max_clinical_stage = r$maxClinicalStage %||%
                                 (drug$maximumClinicalStage %||% NA_character_),
            stringsAsFactors = FALSE)
        if (expand == "mechanism") {
            out <- merge(base, moa, by = NULL)
            out$disease_id   <- collapse(dz$disease_id)
            out$disease_name <- collapse(dz$disease_name)
        } else if (expand == "disease") {
            out <- merge(base, dz, by = NULL)
            out$mechanism_of_action <- collapse(moa$mechanism_of_action)
            out$action_type         <- collapse(moa$action_type)
        } else {  # drug
            out <- base
            out$mechanism_of_action <- collapse(moa$mechanism_of_action)
            out$action_type         <- collapse(moa$action_type)
            out$disease_id          <- collapse(dz$disease_id)
            out$disease_name        <- collapse(dz$disease_name)
        }
        out[, .dtiDrugCols]
    })
    df <- do.call(rbind, perRow)
    rownames(df) <- NULL
    df
}

#' Retrieve known drugs / clinical candidates for one or more targets from
#' Open Targets (target -> drug)
#'
#' Accepts many genes (symbols and/or Ensembl IDs), resolves symbols in
#' bulk via \code{\link{getOpenTargetsIds}}, then pulls
#' \code{drugAndClinicalCandidates} for up to \code{chunkSize} targets per
#' HTTP request using GraphQL aliasing. One input target can expand into
#' several output rows (a drug may act on the target via more than one
#' mechanism and may be developed against more than one disease); the
#' cartesian expansion is controlled by \code{expand}. Note that
#' agonist / undrugged-direct targets (e.g. FGF21, NLRP3, TFEB, ADIPOR)
#' legitimately return zero rows.
#'
#' The Open Targets GraphQL schema evolves between platform releases;
#' this function targets the v4 schema in which rows are of type
#' \code{ClinicalTargetFromTarget} (fields \code{maxClinicalStage},
#' \code{drug}, \code{diseases}) and mechanism data lives under
#' \code{drug.mechanismsOfAction.rows}. If Open Targets changes these
#' names, update \code{\link{.dtiDrugSelection}}.
#'
#' @param genes character vector of gene symbols and/or Ensembl gene IDs.
#' @param expand character(1); \code{"mechanism"} (default) emits one row
#'   per drug x mechanism-of-action and collapses diseases into a single
#'   semicolon-delimited string; \code{"disease"} emits one row per
#'   drug x disease and collapses mechanisms; \code{"drug"} emits one row
#'   per drug with both mechanisms and diseases collapsed.
#' @param chunkSize integer(1) targets per HTTP request (default 25; kept
#'   modest so per-request payloads and server load stay reasonable).
#' @param pause numeric(1) seconds to sleep between requests.
#' @param verbose logical(1) progress messages.
#' @return A \code{data.frame} with columns \code{ensembl_id},
#'   \code{approved_symbol}, \code{drug_id} (ChEMBL), \code{drug_name},
#'   \code{drug_type}, \code{max_clinical_stage}, \code{mechanism_of_action},
#'   \code{action_type}, \code{disease_id}, \code{disease_name}. Genes with
#'   no known drugs simply contribute no rows. Empty \code{data.frame} if
#'   nothing / offline.
#' @examples
#' \donttest{
#'   drugs <- getOpenTargetsDrugs(c("FGFR1", "KLB"))
#'   head(drugs[, c("drug_name", "mechanism_of_action", "max_clinical_stage")])
#' }
#' @seealso \code{\link{getOpenTargetsIds}}, \code{\link{getOpenTargetsTargets}},
#'   \code{\link{getOpenTargetsDrugTarget}}
#' @export
getOpenTargetsDrugs <- function(genes, expand = c("mechanism", "disease", "drug"),
                                chunkSize = 25L, pause = 0.1, verbose = FALSE) {
    expand <- match.arg(expand)
    stopifnot(is.character(genes), length(genes) >= 1L)
    url <- .dtiEndpoints()$opentargets

    ## Resolve any symbols to Ensembl IDs (Ensembl IDs pass through).
    isEnsg <- grepl("^ENSG[0-9]+$", genes)
    ensg <- genes
    if (any(!isEnsg)) {
        mapped <- getOpenTargetsIds(genes[!isEnsg], chunkSize = 100L,
                                    pause = pause, verbose = verbose)
        ensg[!isEnsg] <- mapped
    }
    keep <- !is.na(ensg) & nzchar(ensg)
    ensg <- unique(ensg[keep])
    if (length(ensg) == 0L) return(.dtiEmptyDrugs())

    acc <- vector("list", length(ensg)); k <- 0L
    chunks <- .dtiChunk(ensg, chunkSize)
    for (i in seq_along(chunks)) {
        ch <- chunks[[i]]
        if (verbose) message("getOpenTargetsDrugs chunk ", i, "/", length(chunks),
                             " (", length(ch), " targets)")
        varDefs <- paste(sprintf("$e%d: String!", seq_along(ch)), collapse = ", ")
        aliases <- paste(sprintf("t%d: target(ensemblId: $e%d) { %s }",
                                 seq_along(ch), seq_along(ch),
                                 .dtiDrugSelection()), collapse = "\n")
        query <- sprintf("query drugsBatch(%s) {\n%s\n}", varDefs, aliases)
        vars <- stats::setNames(as.list(ch), sprintf("e%d", seq_along(ch)))
        data <- .dtiGraphQL(url, query, variables = vars)
        if (!is.null(data)) for (j in seq_along(ch)) {
            df <- .dtiParseTargetDrugs(data[[sprintf("t%d", j)]], ch[j], expand)
            if (nrow(df) > 0L) { k <- k + 1L; acc[[k]] <- df }
        }
        if (pause > 0 && i < length(chunks)) Sys.sleep(pause)
    }
    if (k == 0L) return(.dtiEmptyDrugs())
    out <- do.call(rbind, acc[seq_len(k)])
    rownames(out) <- NULL
    out
}

#' Resolve many drug names to ChEMBL IDs in batched GraphQL requests
#'
#' Vectorised, batched the same way \code{\link{getOpenTargetsIds}}
#' batches target-symbol resolution. Open Targets drug IDs *are* ChEMBL
#' IDs.
#'
#' @param names character vector of drug names (ChEMBL IDs pass through
#'   unresolved - checked via a plain \code{"^CHEMBL[0-9]+$"} regex).
#' @param chunkSize integer(1) names per HTTP request (default 100).
#' @param pause numeric(1) seconds to sleep between requests (politeness).
#' @param verbose logical(1) progress messages.
#' @return A named character vector (names = input drug names) of ChEMBL
#'   IDs; unresolved names are \code{NA}. Order matches the input.
#' @examples
#' \donttest{
#'   getOpenTargetsDrugIds(c("aspirin", "imatinib", "NOT_A_DRUG"))
#' }
#' @seealso \code{\link{getOpenTargetsTargets}}
#' @export
getOpenTargetsDrugIds <- function(names, chunkSize = 100L, pause = 0.1,
                                  verbose = FALSE) {
    stopifnot(is.character(names), length(names) >= 1L)
    url <- .dtiEndpoints()$opentargets
    uniq <- unique(names)
    chunks <- .dtiChunk(uniq, chunkSize)
    resolved <- character(0)
    for (i in seq_along(chunks)) {
        ch <- chunks[[i]]
        if (verbose) message("getOpenTargetsDrugIds chunk ", i, "/", length(chunks),
                             " (", length(ch), " names)")
        varDefs <- paste(sprintf("$q%d: String!", seq_along(ch)), collapse = ", ")
        aliases <- paste(sprintf(
            "a%d: search(queryString: $q%d, entityNames: [\"drug\"], page: {index: 0, size: 3}) { hits { id name } }",
            seq_along(ch), seq_along(ch)), collapse = "\n")
        query <- sprintf("query resolveDrugIds(%s) {\n%s\n}", varDefs, aliases)
        vars <- stats::setNames(as.list(ch), sprintf("q%d", seq_along(ch)))
        data <- .dtiGraphQL(url, query, variables = vars)
        for (j in seq_along(ch)) {
            hits <- data[[sprintf("a%d", j)]]$hits %||% list()
            id <- NA_character_
            if (length(hits) > 0L) {
                exact <- Filter(function(h) {
                    nm <- h$name %||% NA_character_
                    !is.na(nm) && toupper(nm) == toupper(ch[j])
                }, hits)
                pick <- if (length(exact) > 0L) exact[[1]] else hits[[1]]
                id <- pick$id %||% NA_character_
            }
            resolved[ch[j]] <- id
        }
        if (pause > 0 && i < length(chunks)) Sys.sleep(pause)
    }
    resolved[names]  # re-expand to input order/length
}

#' Column order shared by the target accessor
#' @keywords internal
.dtiTargetCols <- c("chembl_id", "drug_name", "drug_type", "max_clinical_stage",
                    "mechanism_of_action", "action_type", "moa_target_name",
                    "target_id", "approved_symbol")

#' Empty target data.frame with the canonical columns
#' @keywords internal
.dtiEmptyTargets <- function() {
    as.data.frame(stats::setNames(
        replicate(length(.dtiTargetCols), character(0), simplify = FALSE),
        .dtiTargetCols), stringsAsFactors = FALSE)
}

#' The reusable GraphQL selection for one drug's targets (no outer braces)
#' @keywords internal
.dtiTargetSelection <- function() {
    "id name drugType maximumClinicalStage
     mechanismsOfAction {
       rows {
         mechanismOfAction actionType targetName
         targets { id approvedSymbol }
       }
     }"
}

#' Parse one drug object into a tidy data.frame, applying the requested
#' row expansion.
#'
#' @param drg list; the GraphQL \code{drug} object (may be NULL).
#' @param chemblId character(1); the ChEMBL ID this object was queried with.
#' @param expand character(1); one of "target", "mechanism".
#' @return data.frame with columns \code{.dtiTargetCols}; empty if no rows.
#' @keywords internal
.dtiParseDrugTargets <- function(drg, chemblId, expand) {
    drg  <- drg %||% list()
    rows <- drg$mechanismsOfAction$rows %||% list()
    if (length(rows) == 0L) return(.dtiEmptyTargets())
    nm    <- drg$name %||% NA_character_
    dtype <- drg$drugType %||% NA_character_
    mcs   <- drg$maximumClinicalStage %||% NA_character_
    collapse <- function(x) paste(unique(stats::na.omit(x)), collapse = "; ")
    perRow <- lapply(rows, function(r) {
        tgts <- r$targets %||% list()
        tg <- if (length(tgts) == 0L)
            data.frame(target_id = NA_character_,
                       approved_symbol = NA_character_, stringsAsFactors = FALSE)
        else do.call(rbind, lapply(tgts, function(t) data.frame(
            target_id       = t$id %||% NA_character_,
            approved_symbol = t$approvedSymbol %||% NA_character_,
            stringsAsFactors = FALSE)))
        tg <- unique(tg)
        base <- data.frame(
            chembl_id           = chemblId,
            drug_name            = nm,
            drug_type            = dtype,
            max_clinical_stage   = mcs,
            mechanism_of_action  = r$mechanismOfAction %||% NA_character_,
            action_type          = r$actionType        %||% NA_character_,
            moa_target_name      = r$targetName        %||% NA_character_,
            stringsAsFactors = FALSE)
        if (expand == "target") {
            out <- merge(base, tg, by = NULL)
        } else {  # mechanism
            out <- base
            out$target_id       <- collapse(tg$target_id)
            out$approved_symbol <- collapse(tg$approved_symbol)
        }
        out[, .dtiTargetCols]
    })
    df <- do.call(rbind, perRow)
    rownames(df) <- NULL
    df
}

#' Retrieve targets for one or more drugs from Open Targets (drug -> target)
#'
#' Accepts many drugs (names and/or ChEMBL IDs), resolves names in bulk
#' via \code{\link{getOpenTargetsDrugIds}}, then pulls
#' \code{mechanismsOfAction} for up to \code{chunkSize} drugs per HTTP
#' request using GraphQL aliasing. This is the drug -> target counterpart
#' of \code{\link{getOpenTargetsDrugs}}. A single mechanism-of-action row
#' can list more than one target (e.g. aspirin's "Cyclooxygenase
#' inhibitor" mechanism lists both PTGS1 and PTGS2); the cartesian
#' expansion is controlled by \code{expand}.
#'
#' @param drugs character vector of drug names and/or ChEMBL IDs.
#' @param expand character(1); \code{"target"} (default) emits one row per
#'   mechanism-of-action x target; \code{"mechanism"} emits one row per
#'   mechanism-of-action and collapses targets into a single
#'   semicolon-delimited string.
#' @param chunkSize integer(1) drugs per HTTP request (default 25).
#' @param pause numeric(1) seconds to sleep between requests.
#' @param verbose logical(1) progress messages.
#' @return A \code{data.frame} with columns \code{chembl_id},
#'   \code{drug_name}, \code{drug_type}, \code{max_clinical_stage},
#'   \code{mechanism_of_action}, \code{action_type}, \code{moa_target_name},
#'   \code{target_id} (Ensembl), \code{approved_symbol}. Drugs with no
#'   known targets simply contribute no rows. Empty \code{data.frame} if
#'   nothing / offline.
#' @examples
#' \donttest{
#'   tgts <- getOpenTargetsTargets("aspirin")
#'   tgts[, c("approved_symbol", "mechanism_of_action")]
#' }
#' @seealso \code{\link{getOpenTargetsDrugIds}}, \code{\link{getOpenTargetsDrugs}},
#'   \code{\link{getOpenTargetsDrugTarget}}
#' @export
getOpenTargetsTargets <- function(drugs, expand = c("target", "mechanism"),
                                  chunkSize = 25L, pause = 0.1, verbose = FALSE) {
    expand <- match.arg(expand)
    stopifnot(is.character(drugs), length(drugs) >= 1L)
    url <- .dtiEndpoints()$opentargets

    ## Resolve any names to ChEMBL IDs (ChEMBL IDs pass through).
    isChembl <- grepl("^CHEMBL[0-9]+$", drugs)
    chemblId <- drugs
    if (any(!isChembl)) {
        mapped <- getOpenTargetsDrugIds(drugs[!isChembl], chunkSize = 100L,
                                        pause = pause, verbose = verbose)
        chemblId[!isChembl] <- mapped
    }
    keep <- !is.na(chemblId) & nzchar(chemblId)
    chemblId <- unique(chemblId[keep])
    if (length(chemblId) == 0L) return(.dtiEmptyTargets())

    acc <- vector("list", length(chemblId)); k <- 0L
    chunks <- .dtiChunk(chemblId, chunkSize)
    for (i in seq_along(chunks)) {
        ch <- chunks[[i]]
        if (verbose) message("getOpenTargetsTargets chunk ", i, "/", length(chunks),
                             " (", length(ch), " drugs)")
        varDefs <- paste(sprintf("$c%d: String!", seq_along(ch)), collapse = ", ")
        aliases <- paste(sprintf("d%d: drug(chemblId: $c%d) { %s }",
                                 seq_along(ch), seq_along(ch),
                                 .dtiTargetSelection()), collapse = "\n")
        query <- sprintf("query targetsBatch(%s) {\n%s\n}", varDefs, aliases)
        vars <- stats::setNames(as.list(ch), sprintf("c%d", seq_along(ch)))
        data <- .dtiGraphQL(url, query, variables = vars)
        if (!is.null(data)) for (j in seq_along(ch)) {
            df <- .dtiParseDrugTargets(data[[sprintf("d%d", j)]], ch[j], expand)
            if (nrow(df) > 0L) { k <- k + 1L; acc[[k]] <- df }
        }
        if (pause > 0 && i < length(chunks)) Sys.sleep(pause)
    }
    if (k == 0L) return(.dtiEmptyTargets())
    out <- do.call(rbind, acc[seq_len(k)])
    rownames(out) <- NULL
    out
}

#' Query Open Targets drug/target data via the uniform queryBy interface
#'
#' Thin \code{queryBy}-dispatching wrapper over
#' \code{\link{getOpenTargetsDrugs}} (target -> drug) and
#' \code{\link{getOpenTargetsTargets}} (drug -> target), matching the
#' \code{queryBy = list(molType, idType, ids)} interface used by
#' \code{\link{drugTargetAnnot}}, \code{\link{getChemblDrugTarget}},
#' \code{\link{getPubchemDrugTarget}} and \code{\link{getDgidbDrugTarget}}.
#' Unlike \code{\link{getDgidbDrugTarget}} (which joins on a
#' case-normalized display name), this joins the \code{QueryIDs} column on
#' the stable resolved ID (Ensembl gene ID / ChEMBL ID) that Open Targets
#' already returns verbatim in every output row, so no casing ambiguity
#' arises; the tradeoff is one extra batched ID-resolution call inside the
#' wrapper (see the file header note).
#'
#' @param queryBy named list with components \code{molType}, \code{idType}
#'   and \code{ids}. \code{molType = "gene"} with \code{idType = "symbol"}
#'   queries target -> drug (\code{ids} = gene symbols and/or Ensembl gene
#'   IDs). \code{molType = "cmp"} with \code{idType = "name"} queries
#'   drug -> target (\code{ids} = drug names and/or ChEMBL IDs). Matches
#'   \code{\link{getPubchemDrugTarget}}/\code{\link{getDgidbDrugTarget}}'s
#'   vocabulary for the same concepts.
#' @param ... additional arguments passed through to
#'   \code{\link{getOpenTargetsDrugs}} / \code{\link{getOpenTargetsTargets}}
#'   (e.g. \code{expand}, \code{chunkSize}, \code{verbose}).
#' @return A \code{data.frame} in the same column shape as
#'   \code{\link{getOpenTargetsDrugs}} / \code{\link{getOpenTargetsTargets}},
#'   plus a leading \code{QueryIDs} column echoing the original
#'   \code{queryBy$ids} token each row resolved from. Query IDs that
#'   returned no rows still appear as a single row with all other fields
#'   \code{NA}.
#' @examples
#' \donttest{
#'   ## target -> drug: FGFR1, KLB
#'   getOpenTargetsDrugTarget(list(molType = "gene", idType = "symbol",
#'                                 ids = c("FGFR1", "KLB")))
#'   ## drug -> target: aspirin
#'   getOpenTargetsDrugTarget(list(molType = "cmp", idType = "name",
#'                                 ids = "aspirin"))
#' }
#' @seealso \code{\link{getOpenTargetsDrugs}}, \code{\link{getOpenTargetsTargets}},
#'   \code{\link{getChemblDrugTarget}}, \code{\link{getPubchemDrugTarget}},
#'   \code{\link{getDgidbDrugTarget}}
#' @export
getOpenTargetsDrugTarget <- function(queryBy = list(molType = NULL, idType = NULL,
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
            "getOpenTargetsDrugTarget() currently supports only ",
            "queryBy=list(molType=\"gene\", idType=\"symbol\", ids=...) ",
            "or queryBy=list(molType=\"cmp\", idType=\"name\", ids=...). ",
            "Other identifier types require translating to a gene symbol ",
            "/ Ensembl ID or a drug name / ChEMBL ID first."
        )
    }

    ids <- queryBy$ids
    if (isGene) {
        out <- getOpenTargetsDrugs(ids, ...)
        resolved <- ids
        isNative <- grepl("^ENSG[0-9]+$", ids)
        if (any(!isNative))
            resolved[!isNative] <- getOpenTargetsIds(ids[!isNative])
        resolvedCol <- "ensembl_id"
    } else {
        out <- getOpenTargetsTargets(ids, ...)
        resolved <- ids
        isNative <- grepl("^CHEMBL[0-9]+$", ids)
        if (any(!isNative))
            resolved[!isNative] <- getOpenTargetsDrugIds(ids[!isNative])
        resolvedCol <- "chembl_id"
    }

    ## Join on the stable resolved ID, not a display name (see file header
    ## note) - `resolved` and `ids` are parallel/same-length, so the first
    ## query token whose resolved ID matches a given output row wins if
    ## more than one input token resolves to the same ID (a rare synonym
    ## case, not fixed here - same tradeoff already accepted for the
    ## ChEMBL/PubChem/DGIdb wrappers' analogous edge cases).
    out$QueryIDs <- ids[match(out[[resolvedCol]], resolved)]

    unmatched <- ids[is.na(resolved) | !(resolved %in% out[[resolvedCol]])]
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
