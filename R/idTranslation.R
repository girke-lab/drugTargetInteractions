## =====================================================================
##  idTranslation.R
##  Generic ID-translation layer for the drugTargetInteractions
##  Bioconductor package: protein/gene ID mapping (UniProt REST) and
##  paralog/ortholog lookup (Ensembl REST homology).
##
##  Both go direct against the underlying public REST APIs rather than
##  through the R wrapper packages this package already depends on
##  (UniProt.ws, biomaRt), decided 2026-07-17:
##   1. UniProt.ws is a thin wrapper around this exact REST ID-mapping
##      service and has a multi-year history of lagging UniProt's own
##      API changes (see getUniprotIDs()'s existing breakage - keytype/
##      column renames, and a real correctness bug where mapUniProt()'s
##      taxId organism filter was silently non-functional for roughly a
##      year before UniProt.ws 2.49.1). The raw REST API's taxId filter
##      was confirmed live to work correctly.
##   2. biomaRt/Ensembl BioMart has its own long-documented reliability
##      problems (timeouts, intermittent outages) independent of
##      UniProt.ws, including reports specifically about paralog queries
##      timing out. The Ensembl REST /homology/ endpoint serves the same
##      underlying Compara data with a lean, rate-documented (15 req/s)
##      interface, and additionally exposes true cross-species orthologs
##      (not just within-human paralogs, unlike the existing
##      getParalogs()).
##
##  getUniprotIDs() (UniProt.ws) and getParalogs() (biomaRt) are left in
##  place, unmodified, as historical/backup options - not removed, not
##  extended, not the default recommended path going forward. Worth a
##  "legacy/alternative approach" mention in a future vignette supplement.
##
##  Ensembl's GOC/WGA-based "high confidence ortholog" scoring
##  (orthology_confidence/goc_score/wga_coverage) is real - it combines
##  sequence identity with gene-order conservation and whole-genome-
##  alignment coverage - but is only reachable via BioMart, not this REST
##  endpoint (confirmed by live testing and by grepping the ensembl-rest
##  server source). Deliberately not chased here since it would
##  reintroduce the exact BioMart reliability problem this module avoids,
##  for orthologs only (paralogs never had it either way - GOC/WGA are
##  inherently cross-genome synteny concepts). Ensembl's homology `type`
##  classification (ortholog_one2one/one2many/many2many, other_paralog,
##  within_species_paralog, ...) is already a gene-tree-reconciliation
##  call, not raw sequence-identity thresholding, and is the recommended
##  "best hit" filter here: prefer `type == "ortholog_one2one"` (or the
##  lowest-fan-out paralog `type`), then rank ties by `perc_id`.
## =====================================================================


## ---------------------------------------------------------------------
## UniProt REST ID Mapping (replaces UniProt.ws for translation purposes)
## ---------------------------------------------------------------------
## Job-based: POST /idmapping/run -> poll /idmapping/status/{jobId} ->
## paginated GET /idmapping/results/{jobId} (follows the `Link: rel=
## "next"` header, confirmed live to appear once a job's results exceed
## one page). `from`/`to` use UniProt's own database-name vocabulary
## directly (e.g. "Gene_Name", "UniProtKB-Swiss-Prot", "Ensembl",
## "UniProtKB_AC-ID") rather than a package-specific abstraction, since
## that vocabulary is itself UniProt's stable, versioned public contract
## - see https://www.uniprot.org/help/id_mapping for the full list.

#' Submit a UniProt ID mapping job
#' @keywords internal
#' @noRd
.dtiUniprotSubmitJob <- function(ids, from, to, taxId = NULL) {
    url <- paste0(.dtiEndpoints()$uniprot, "/idmapping/run")
    body <- list(from = from, to = to, ids = paste(unique(ids), collapse = ","))
    if (!is.null(taxId)) body$taxId <- as.character(taxId)
    res <- .dtiApiPOSTform(url, body, hardStop = TRUE)
    res$jobId %||% stop("UniProt ID mapping job submission did not return a jobId.")
}

#' Poll a UniProt ID mapping job until finished
#'
#' UniProt's ID mapping is asynchronous even for tiny jobs (confirmed
#' live: a single-ID job is still "RUNNING" immediately after
#' submission), so polling is always needed - there is no synchronous
#' fast path.
#'
#' The status endpoint responds \code{303 See Other} once the job is
#' done, with the redirect \code{Location} pointing at a results view
#' and the literal \code{\{"jobStatus":"FINISHED"\}} body attached to
#' the 303 response itself. httr2 follows redirects by default like
#' most HTTP clients, so a plain GET here would silently land on the
#' results page instead of the status payload - confirmed live to hang
#' forever otherwise, since \code{jobStatus} is then never present to
#' detect. Redirect-following is disabled explicitly for this request
#' only (not changed package-wide in \code{.dtiApiGET}).
#' @keywords internal
#' @noRd
.dtiUniprotPollJob <- function(jobId, pollInterval = 1, maxWait = 120) {
    url <- paste0(.dtiEndpoints()$uniprot, "/idmapping/status/", jobId)
    waited <- 0
    repeat {
        res <- tryCatch({
            req <- httr2::request(url)
            req <- httr2::req_headers(req, Accept = "application/json")
            req <- httr2::req_user_agent(
                req, "drugTargetInteractions R package (Bioconductor)")
            req <- httr2::req_timeout(req, 60L)
            req <- httr2::req_options(req, followlocation = 0L)
            req <- httr2::req_throttle(req, rate = 5, fill_time_s = 1)
            req <- httr2::req_retry(
                req, max_tries = 3L,
                is_transient = function(resp)
                    httr2::resp_status(resp) %in% c(429L, 500L, 502L, 503L, 504L))
            resp <- httr2::req_perform(req)
            httr2::resp_body_json(resp, simplifyVector = FALSE)
        }, error = function(e) stop(e))
        status <- res$jobStatus %||% NA_character_
        if (identical(status, "FINISHED")) return(invisible(TRUE))
        if (status %in% c("ERROR", "FAILED"))
            stop("UniProt ID mapping job '", jobId, "' failed.")
        if (waited >= maxWait)
            stop("UniProt ID mapping job '", jobId, "' timed out after ",
                 maxWait, "s.")
        Sys.sleep(pollInterval)
        waited <- waited + pollInterval
    }
}

#' Fetch all pages of a finished UniProt ID mapping job's results
#'
#' Pages via the response's \code{Link: rel="next"} header. Treats a 404
#' as transient-and-retryable on top of the usual 429/5xx set: live
#' testing found a short race window immediately after a job's status
#' turns \code{"FINISHED"} where the results endpoint can still 404 once
#' before succeeding.
#' @keywords internal
#' @noRd
.dtiUniprotFetchResults <- function(jobId, pageSize = 500L) {
    url <- paste0(.dtiEndpoints()$uniprot, "/idmapping/results/", jobId)
    acc <- list()
    nextUrl <- url
    nextQuery <- list(format = "json", size = pageSize)
    repeat {
        req <- httr2::request(nextUrl)
        req <- httr2::req_headers(req, Accept = "application/json")
        req <- httr2::req_user_agent(
            req, "drugTargetInteractions R package (Bioconductor)")
        req <- httr2::req_timeout(req, 60L)
        if (!is.null(nextQuery)) req <- httr2::req_url_query(req, !!!nextQuery)
        req <- httr2::req_throttle(req, rate = 5, fill_time_s = 1)
        req <- httr2::req_retry(
            req, max_tries = 3L,
            is_transient = function(resp)
                httr2::resp_status(resp) %in% c(404L, 429L, 500L, 502L, 503L, 504L))
        resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
        if (is.null(resp)) break
        parsed <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        acc <- c(acc, parsed$results %||% list())
        nxt <- tryCatch(httr2::resp_link_url(resp, "next"), error = function(e) NULL)
        if (is.null(nxt)) break
        nextUrl <- nxt
        nextQuery <- NULL  ## already embedded in the Link URL
    }
    acc
}

#' Tidy a UniProt ID mapping job's raw results into a data.frame
#' @keywords internal
#' @noRd
.dtiUniprotParseResults <- function(results) {
    cols <- c("From", "To")
    empty <- as.data.frame(stats::setNames(
        replicate(length(cols), character(0), simplify = FALSE), cols),
        stringsAsFactors = FALSE)
    if (length(results) == 0L) return(empty)
    df <- do.call(rbind, lapply(results, function(r) data.frame(
        From = r$from %||% NA_character_,
        To   = if (is.list(r$to)) (r$to$primaryAccession %||% r$to$id %||% NA_character_)
               else (r$to %||% NA_character_),
        stringsAsFactors = FALSE)))
    rownames(df) <- NULL
    df
}

#' Map identifiers between databases via the UniProt REST ID Mapping API
#'
#' Direct replacement for \code{\link{getUniprotIDs}} (UniProt.ws), which
#' is left in place unmodified as a historical/backup option but is not
#' the recommended path going forward - see the \code{idTranslation.R}
#' file header for the full rationale (reliability and a real
#' organism-filtering bug in UniProt.ws that this direct implementation
#' avoids by construction, since it calls the same REST service
#' UniProt.ws itself wraps).
#'
#' Submits an asynchronous ID mapping job, polls until it finishes, and
#' returns every \code{(From, To)} pair found, paginating through all
#' result pages automatically.
#'
#' @param ids character vector of source-database identifiers (e.g. gene
#'   symbols, UniProt accessions, Ensembl gene IDs). Deduplicated
#'   internally before submission.
#' @param from character(1) source database name, in UniProt's own
#'   vocabulary (e.g. \code{"Gene_Name"}, \code{"Ensembl"},
#'   \code{"UniProtKB_AC-ID"}). See
#'   \url{https://www.uniprot.org/help/id_mapping} for the full list.
#' @param to character(1) target database name (e.g.
#'   \code{"UniProtKB-Swiss-Prot"} for canonical reviewed entries only,
#'   \code{"UniProtKB"} for reviewed+unreviewed, \code{"Ensembl"},
#'   \code{"Gene_Name"}).
#' @param taxId integer(1) or \code{NULL}; NCBI taxonomy ID to restrict
#'   matches to one organism (default \code{9606L} = human), or
#'   \code{NULL} to span all organisms matching \code{ids}. Applied
#'   server-side at job-submission time (confirmed live to filter
#'   correctly, unlike the historical UniProt.ws bug noted above).
#' @param pollInterval numeric(1) seconds between job-status polls.
#' @param maxWait numeric(1) seconds to wait for job completion before
#'   erroring.
#' @param verbose logical(1); if TRUE, message the job ID and size on
#'   submission.
#' @return A \code{data.frame} with columns \code{From} and \code{To}.
#'   IDs that don't map contribute no row (matching the underlying API's
#'   behavior - there is no per-ID failure signal to NA-pad against, only
#'   an aggregate list of unmapped IDs the API does not currently surface
#'   in the paginated JSON results). Note: \code{to = "Ensembl"} results
#'   carry a version suffix (e.g. \code{"ENSG00000077782.24"}), unlike
#'   the bare IDs returned elsewhere in this package - strip
#'   \code{sub("\\\\..*$", "", x)} if a bare ID is needed.
#' @examples
#' \donttest{
#'   getUniprotMapping(c("FGFR1", "KLB"), from = "Gene_Name",
#'                     to = "UniProtKB-Swiss-Prot")
#'   getUniprotMapping("P11362", from = "UniProtKB_AC-ID", to = "Ensembl")
#' }
#' @seealso \code{\link{getUniprotIDs}}
#' @export
getUniprotMapping <- function(ids, from, to, taxId = 9606L, pollInterval = 1,
                              maxWait = 120, verbose = FALSE) {
    stopifnot(is.character(ids), length(ids) >= 1L,
             is.character(from), length(from) == 1L,
             is.character(to), length(to) == 1L)
    jobId <- .dtiUniprotSubmitJob(ids, from = from, to = to, taxId = taxId)
    if (verbose)
        message("UniProt ID mapping job '", jobId, "' submitted (", from,
                " -> ", to, ", ", length(unique(ids)), " IDs)")
    .dtiUniprotPollJob(jobId, pollInterval = pollInterval, maxWait = maxWait)
    .dtiUniprotParseResults(.dtiUniprotFetchResults(jobId))
}


## ---------------------------------------------------------------------
## Ensembl REST homology (paralogs + orthologs, replaces biomaRt for
## translation purposes)
## ---------------------------------------------------------------------
## The endpoint's default response embeds a full pairwise protein
## alignment (`cigar_line` + `align_seq`) per homolog, which is fine for
## an isolated gene but pathological for large gene families - live-
## tested 2026-07-17: FGFR1 (53 paralogs, part of the RTK superfamily)
## took 72s uncapped vs. 6s with `cigar_line=0` (drops the CIGAR string,
## keeps perc_id/perc_pos) vs. 1.4s with `format=condensed` (drops
## perc_id too, gene IDs + relationship type only). `cigar_line=0` is
## the default query shape here; pass `condensed = TRUE` for the faster,
## ranking-free existence-check path.

#' Column order for tidy Ensembl homology rows
#' @keywords internal
#' @noRd
.dtiEnsemblHomologyCols <- c("query_gene", "homolog_id", "homolog_protein_id",
                             "homolog_species", "type", "taxonomy_level",
                             "perc_id", "perc_pos")

#' Empty Ensembl homology data.frame with the canonical columns
#' @keywords internal
#' @noRd
.dtiEmptyEnsemblHomology <- function() {
    as.data.frame(stats::setNames(
        replicate(length(.dtiEnsemblHomologyCols), character(0), simplify = FALSE),
        .dtiEnsemblHomologyCols), stringsAsFactors = FALSE)
}

#' Fetch + parse homology records for one gene
#'
#' Shared by \code{\link{getEnsemblParalogs}} and
#' \code{\link{getEnsemblOrthologs}}, which only differ in \code{type}.
#' @keywords internal
#' @noRd
.dtiEnsemblHomologyFetch <- function(gene, species = "human",
                                     type = c("paralogues", "orthologues"),
                                     targetSpecies = NULL, condensed = FALSE,
                                     verbose = FALSE) {
    type <- match.arg(type)
    url <- paste0(.dtiEndpoints()$ensembl, "/homology/symbol/", species, "/",
                 utils::URLencode(gene, reserved = TRUE))
    q <- list(type = type)
    if (!is.null(targetSpecies)) q$target_species <- targetSpecies
    if (condensed) q$format <- "condensed" else q$cigar_line <- 0L
    if (verbose) message("Ensembl homology (", type, ") for ", gene)
    ## Ensembl REST latency varies a lot with server load (live-observed:
    ## an isolated-gene query that took ~1.3s in one session took ~5.8s
    ## in another), and large gene families compound that - a generous
    ## timeout here is a deliberate tolerance for that variability, not a
    ## sign the default 60s in .dtiApiGET() is normally too short.
    res <- .dtiApiGET(url, query = q, timeout = 120L)
    homs <- res$data[[1]]$homologies %||% list()
    if (length(homs) == 0L) return(.dtiEmptyEnsemblHomology())
    ## `format=condensed` puts id/protein_id/species at the top level of
    ## each homology record; the default (cigar_line=0) format nests
    ## them under `target` alongside perc_id/perc_pos instead. Fall back
    ## to `h` itself when `target` is absent so both shapes parse.
    out <- do.call(rbind, lapply(homs, function(h) {
        tgt <- h$target %||% h
        data.frame(
            query_gene         = gene,
            homolog_id         = tgt$id %||% NA_character_,
            homolog_protein_id = tgt$protein_id %||% NA_character_,
            homolog_species    = tgt$species %||% NA_character_,
            type               = h$type %||% NA_character_,
            taxonomy_level     = h$taxonomy_level %||% NA_character_,
            perc_id            = as.numeric(tgt$perc_id %||% NA),
            perc_pos           = as.numeric(tgt$perc_pos %||% NA),
            stringsAsFactors = FALSE)
    }))
    rownames(out) <- NULL
    out
}

#' Retrieve within-species paralogs for one or more genes from Ensembl
#'
#' Direct replacement for \code{\link{getParalogs}} (biomaRt), which is
#' left in place unmodified as a historical/backup option but is not the
#' recommended path going forward - see the \code{idTranslation.R} file
#' header for the full rationale.
#'
#' \code{type} (e.g. \code{"other_paralog"}, \code{"within_species_
#' paralog"}, \code{"gene_split"}) and \code{taxonomy_level} come from
#' Ensembl Compara's gene-tree reconciliation, not raw sequence-identity
#' thresholding - already a more principled "relatedness" signal than
#' \code{perc_id} alone. \code{perc_id}/\code{perc_pos} are useful for
#' ranking within a \code{type} tier.
#'
#' @param genes character vector of gene symbols.
#' @param species character(1) Ensembl species name or alias for
#'   \code{genes} (default \code{"human"}).
#' @param condensed logical(1); if \code{TRUE}, use the faster
#'   \code{format=condensed} response shape (gene/protein IDs + \code{type}
#'   only, no \code{perc_id}/\code{perc_pos}) - much faster for
#'   large gene families (see file header), useful when only an
#'   existence check or the relationship \code{type} is needed.
#' @param verbose logical(1); if TRUE, message progress per gene.
#' @return A \code{data.frame} with columns \code{query_gene},
#'   \code{homolog_id}, \code{homolog_protein_id}, \code{homolog_species}
#'   (always the same as \code{species} here), \code{type},
#'   \code{taxonomy_level}, \code{perc_id}, \code{perc_pos} (the latter
#'   two \code{NA} when \code{condensed = TRUE}). Genes with no paralogs
#'   simply contribute no rows.
#' @examples
#' \donttest{
#'   getEnsemblParalogs("FGFR1")
#' }
#' @seealso \code{\link{getEnsemblOrthologs}}, \code{\link{getParalogs}}
#' @export
getEnsemblParalogs <- function(genes, species = "human", condensed = FALSE,
                               verbose = FALSE) {
    stopifnot(is.character(genes), length(genes) >= 1L)
    rows <- lapply(genes, .dtiEnsemblHomologyFetch, species = species,
                   type = "paralogues", condensed = condensed, verbose = verbose)
    rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
    if (length(rows) == 0L) return(.dtiEmptyEnsemblHomology())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' Retrieve cross-species orthologs for one or more genes from Ensembl
#'
#' Unlike \code{\link{getParalogs}} (biomaRt, human paralogs only), this
#' also covers true cross-species orthologs - see
#' \code{\link{getEnsemblParalogs}} for the shared rationale and the
#' \code{type}/\code{taxonomy_level} interpretation.
#'
#' @param genes character vector of gene symbols.
#' @param species character(1) Ensembl species name or alias for
#'   \code{genes} (default \code{"human"}).
#' @param targetSpecies character vector of Ensembl species names/aliases
#'   to restrict orthologs to (e.g. \code{"mouse"}), or \code{NULL}
#'   (default) for all species with a called ortholog.
#' @param condensed logical(1); see \code{\link{getEnsemblParalogs}}.
#' @param verbose logical(1); if TRUE, message progress per gene.
#' @return A \code{data.frame}; see \code{\link{getEnsemblParalogs}} for
#'   columns (\code{homolog_species} varies here). Genes with no
#'   orthologs simply contribute no rows.
#' @examples
#' \donttest{
#'   getEnsemblOrthologs("NLRP3", targetSpecies = "mouse")
#' }
#' @seealso \code{\link{getEnsemblParalogs}}, \code{\link{getParalogs}}
#' @export
getEnsemblOrthologs <- function(genes, species = "human", targetSpecies = NULL,
                                condensed = FALSE, verbose = FALSE) {
    stopifnot(is.character(genes), length(genes) >= 1L)
    rows <- lapply(genes, .dtiEnsemblHomologyFetch, species = species,
                   type = "orthologues", targetSpecies = targetSpecies,
                   condensed = condensed, verbose = verbose)
    rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
    if (length(rows) == 0L) return(.dtiEmptyEnsemblHomology())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
