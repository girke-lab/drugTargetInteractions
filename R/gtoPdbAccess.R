## =====================================================================
##  gtoPdbAccess.R
##  IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) integration for the
##  drugTargetInteractions Bioconductor package.
##
##  GtoPdb publishes its data two ways: a set of public bulk download
##  files under /DATA/, and a REST API under /services/. This module
##  builds from the bulk files. Both the whole interactions table
##  (interactions.tsv, ~25,000 rows) and the small target ID <-> HGNC
##  symbol mapping (GtP_to_HGNC_mapping.tsv, ~3,800 rows) are downloaded
##  once, cached via the package's existing BiocFileCache helpers, and
##  built into a local SQLite - the same build-once-query-many pattern as
##  ttdAccess.R/broadRepurposingHubAccess.R.
##
##  WHY THE BULK FILES AND NOT /services/: as of 2026 GtoPdb requires a
##  registered API key for every /services/ endpoint - an unauthenticated
##  request returns HTTP 401 '{"message":"API key is missing..."}' rather
##  than data (live-confirmed 2026-10-13). This module previously used
##  the bulk /services/interactions endpoint, which meant the package
##  could not be built or checked on any machine without a key,
##  including Bioconductor's own build machines. The /DATA/ bulk files
##  need no key, carry the same release (version comment in their first
##  line), and interactions.tsv actually exposes more columns than the
##  REST payload did - including the PubMed ID of each interaction's
##  reference, where the REST payload carried only GtoPdb's internal
##  reference IDs. The authenticated REST route is still available via
##  downloadGtoPdb(useServices = TRUE, apiKey = ...) for anyone who has
##  a key and prefers it; it is normalized to exactly the same columns,
##  except that it cannot populate pubmedId.
##
##  A real data-modeling gotcha in GtP_to_HGNC_mapping.tsv: its "IUPHAR
##  ID" column mixes TWO separate ID namespaces (target/object IDs and
##  ligand IDs) with no type flag - they collide numerically (e.g. ID 684
##  is both the ligand "ADM2" and the target "ASIC1"). The row's "GtP
##  URL" column disambiguates them (.../ObjectDisplayForward?... for
##  targets vs .../LigandDisplayForward?... for ligands) - this module
##  filters to ObjectDisplayForward rows only before joining against
##  interactions$targetId. A target ID legitimately maps to more than one
##  HGNC symbol for multi-subunit complexes (e.g. the F-type ATPase C
##  subunit, TargetID 803, maps to ATP5MC1/ATP5MC2/ATP5MC3) - kept as
##  separate rows, not collapsed, mirroring
##  broadRepurposingHubAccess.R's .brhExplodeTargets() philosophy of one
##  row per (interaction, gene) pair.
##
##  interactions.tsv also carries 572 rows (298 of them human) with an
##  empty "Target ID", "Target" and "Target Gene Symbol" - ligands whose
##  target is not a GtoPdb object record. They cannot form a
##  (target, drug) edge, which is the entire content of this table, so
##  they are dropped at build time rather than kept as rows with a NA
##  target.
##
##  Interactions are filtered to targetSpecies == "Human" - the rest of
##  this package is implicitly human-gene-centric (HGNC anchor, human
##  UniProt accessions), and GtP_to_HGNC_mapping.tsv itself mixes
##  human/mouse/rat gene symbols in its "HGNC Symbol" column for
##  cross-species target entries, so species-filtering the interactions
##  side first avoids any risk of a cross-species symbol collision.
##
##  LICENSE / distribution posture: GtoPdb data is licensed under the
##  Open Data Commons Open Database License (ODbL) with contents under
##  CC-BY-SA 4.0 - clear terms that explicitly permit redistribution
##  (with attribution/share-alike), unlike TTD's ambiguous or the Broad
##  Repurposing Hub's non-commercial-only wording. That is what makes
##  the small inst/extdata/gtopdb_interactions_example.tsv extract
##  shipped for the examples and the vignette possible at all; the full
##  table is still built into the caller's own cache at call time, as
##  for the other local-SQLite sources.
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for GtoPdb's public bulk download files
#' @keywords internal
#' @noRd
.gtoPdbEndpoints <- function() {
    list(
        interactions = "https://www.guidetopharmacology.org/DATA/interactions.tsv",
        hgncMapping  = "https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.tsv"
    )
}

#' GtoPdb's REST interactions endpoint. Requires a registered API key
#' (see this file's header) - only reached via
#' \code{downloadGtoPdb(useServices = TRUE)}.
#' @keywords internal
#' @noRd
.gtoPdbServicesUrl <- function() {
    "https://www.guidetopharmacology.org/services/interactions"
}

#' Cache resource names, kept distinct per route so a flat-file copy and
#' a REST copy of the interactions data can never overwrite each other.
#' @keywords internal
#' @noRd
.gtoPdbCacheNames <- function(useServices = FALSE) {
    c(interactions = if (isTRUE(useServices)) "gtp_interactions.json"
                     else "gtp_interactions.tsv",
      hgncMapping  = "GtP_to_HGNC_mapping.tsv")
}

#' Strip HTML markup (GtoPdb target/ligand names embed \\code{<sub>}/
#' \\code{<sup>}/\\code{<i>}/\\code{<small>} tags, e.g.
#' \\code{"5-HT<sub>1A</sub> receptor"}) down to plain text.
#' @keywords internal
#' @noRd
.gtoPdbStripHtml <- function(x) {
    gsub("<[^>]+>", "", x)
}

#' Extract GtoPdb's own release version from the
#' "# GtoPdb Version: X.Y - published: YYYY-MM-DD" first-line comment
#' that every /DATA/ bulk file carries.
#' @keywords internal
#' @noRd
.gtoPdbVersion <- function(path) {
    hdr <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
    hit <- regmatches(hdr, regexpr("Version:\\s*[0-9.]+", hdr))
    if (length(hit) == 0L) return(NA_character_)
    trimws(sub("Version:\\s*", "", hit))
}

#' Empty strings to NA. GtoPdb's bulk files leave an unknown value blank
#' where the REST payload sent JSON null; this keeps both routes
#' agreeing on NA.
#' @keywords internal
#' @noRd
.gtoPdbBlankToNa <- function(x) {
    x[!nzchar(trimws(x))] <- NA_character_
    x
}

#' Format one affinity bound the way GtoPdb's own REST payload did: a
#' single decimal place. The bulk file stores these as single-precision
#' floats, so a raw value reads 9.399999618530273 where GtoPdb reports
#' "9.4".
#'
#' Blanks are mapped to NA before the numeric conversion: every value in
#' the three affinity columns is either blank or a number (checked across
#' the whole 2026.3.1 release), and as.numeric(NA_character_) is silent
#' where as.numeric("") warns. Should GtoPdb ever put something else
#' there, the resulting coercion warning is worth seeing rather than
#' suppressing.
#' @keywords internal
#' @noRd
.gtoPdbNum1 <- function(x) {
    v <- as.numeric(.gtoPdbBlankToNa(x))
    ifelse(is.na(v), NA_character_, sprintf("%.1f", v))
}

#' Rebuild the REST payload's single \code{affinity} string from the
#' bulk file's separate low/median/high columns: a range
#' ("5.0 - 5.5") where the bounds differ, otherwise the single value
#' ("9.3"). Verified against a cached /services/interactions response:
#' the rule reproduces it exactly on rows the two can be joined
#' unambiguously.
#' @keywords internal
#' @noRd
.gtoPdbAffinity <- function(low, median, high) {
    lo <- .gtoPdbNum1(low)
    md <- .gtoPdbNum1(median)
    hi <- .gtoPdbNum1(high)
    ifelse(!is.na(lo) & !is.na(hi) & lo != hi,
           paste(lo, "-", hi),
           ifelse(!is.na(md), md, ifelse(!is.na(hi), hi, lo)))
}


## ---------------------------------------------------------------------
## Download raw data (cached via the package's existing BiocFileCache
## helpers)
## ---------------------------------------------------------------------

#' Fetch GtoPdb's authenticated REST interactions endpoint
#'
#' Only reached via \code{downloadGtoPdb(useServices = TRUE)}. Kept
#' separate from \code{.downloadFile()} because the key travels in a
#' request header, which BiocFileCache's own URL-based download cannot
#' attach; the response is written to a temp file and then handed to the
#' cache. HTTP status is inspected here rather than left to httr2's
#' default so a rejected or missing key produces a message naming the
#' cause.
#' @keywords internal
#' @noRd
.gtoPdbFetchServices <- function(apiKey, name) {
    if (!is.character(apiKey) || length(apiKey) != 1L || is.na(apiKey) ||
        !nzchar(apiKey))
        stop("GtoPdb's REST web services require a registered API key. ",
             "Supply one via the 'apiKey' argument or the GTOPDB_API_KEY ",
             "environment variable (request a key at ",
             "https://www.guidetopharmacology.org/webServices.jsp), or ",
             "leave 'useServices = FALSE' to build from GtoPdb's public ",
             "bulk download files instead, which need no key.",
             call. = FALSE)

    url <- .gtoPdbServicesUrl()
    tmp <- tempfile(fileext = ".json")
    req <- httr2::request(url)
    req <- httr2::req_headers(req, `GTP-API-Key` = apiKey,
                              Accept = "application/json")
    req <- httr2::req_user_agent(
        req, "drugTargetInteractions R package (Bioconductor)")
    req <- httr2::req_timeout(req, 600L)
    ## Handle the status here instead of letting httr2 throw its generic
    ## HTTP error, so 401/403 can say what to do about it.
    req <- httr2::req_error(req, is_error = function(resp) FALSE)
    resp <- httr2::req_perform(req, path = tmp)
    status <- httr2::resp_status(resp)
    if (status %in% c(401L, 403L)) {
        unlink(tmp)
        stop("GtoPdb rejected the request to '", url, "' with HTTP ", status,
             " (authentication). Check the API key passed via 'apiKey' or ",
             "the GTOPDB_API_KEY environment variable, or leave ",
             "'useServices = FALSE' to build from GtoPdb's public bulk ",
             "download files instead, which need no key.", call. = FALSE)
    }
    if (status >= 400L) {
        unlink(tmp)
        stop("GtoPdb request to '", url, "' failed with HTTP ", status, ".",
             call. = FALSE)
    }

    out <- .dtiCacheLocalFile(name, tmp)
    unlink(tmp)
    out
}

#' Download GtoPdb's bulk interactions and HGNC mapping data
#'
#' Downloads (or reuses previously cached copies of) the two GtoPdb bulk
#' files this package builds from - the whole interactions table and the
#' target-ID-to-HGNC-symbol mapping - via the package's existing
#' \code{.downloadFile()}/BiocFileCache infrastructure, the same
#' mechanism \code{downloadTTD()}/\code{downloadBroadRepurposingHub()}
#' use. Files land in the local BiocFileCache (see \code{.getCache()},
#' \code{rappdirs::user_cache_dir(appname = "drugTargetInteractions")}) -
#' no GtoPdb data is downloaded until this is called.
#'
#' Neither file needs an API key. GtoPdb's REST web services do, as of
#' 2026, for every \code{/services/} endpoint; set
#' \code{useServices = TRUE} to take the interactions data from
#' \code{/services/interactions} instead and supply \code{apiKey}. A
#' missing or rejected key is reported as such.
#'
#' @param rerun logical(1); if \code{TRUE} (default), check for updates
#'   and (re)download as needed; if \code{FALSE}, use whatever is
#'   already in the local cache without checking upstream.
#' @param config list as returned by \code{genConfig()}; unused here but
#'   accepted for consistency with the rest of the package's API.
#' @param useServices logical(1); if \code{TRUE}, take the interactions
#'   data from GtoPdb's authenticated \code{/services/interactions} REST
#'   endpoint rather than the public bulk file. Requires \code{apiKey}.
#' @param apiKey character(1) GtoPdb API key, used only when
#'   \code{useServices = TRUE}. Defaults to the \code{GTOPDB_API_KEY}
#'   environment variable.
#' @return A named list of local file paths: \code{interactions},
#'   \code{hgncMapping}.
#' @examples
#' ## Downloaded files are cached here and reused by later calls:
#' rappdirs::user_cache_dir("drugTargetInteractions")
#'
#' ## A small extract of both files ships with the package, for the
#' ## examples and the vignette (see buildGtoPdbDb()'s 'paths'):
#' basename(c(
#'     system.file("extdata", "gtopdb_interactions_example.tsv",
#'                 package = "drugTargetInteractions"),
#'     system.file("extdata", "gtopdb_hgnc_mapping_example.tsv",
#'                 package = "drugTargetInteractions")))
#'
#' \donttest{
#'   paths <- downloadGtoPdb()
#'   paths
#' }
#' @seealso \code{\link{buildGtoPdbDb}}
#' @export
downloadGtoPdb <- function(rerun = TRUE, config = genConfig(),
                           useServices = FALSE,
                           apiKey = Sys.getenv("GTOPDB_API_KEY")) {
    ep <- .gtoPdbEndpoints()
    fnames <- .gtoPdbCacheNames(useServices)
    paths <- vector("list", length(ep))
    names(paths) <- names(ep)
    for (nm in names(ep)) {
        if (!rerun) {
            paths[[nm]] <- .getCacheFile(fnames[[nm]])
        } else if (isTRUE(useServices) && nm == "interactions") {
            paths[[nm]] <- .gtoPdbFetchServices(apiKey, fnames[[nm]])
        } else {
            paths[[nm]] <- .downloadFile(ep[[nm]], fnames[[nm]])
        }
    }
    paths
}


## ---------------------------------------------------------------------
## Parse the raw data
## ---------------------------------------------------------------------

#' Parse GtP_to_HGNC_mapping.tsv into a target-ID -> HGNC-symbol lookup,
#' filtered to target/object rows only (see file header for why - the
#' raw file's "IUPHAR ID" column also carries unrelated ligand IDs in
#' the same numeric space).
#' @keywords internal
#' @noRd
.gtoPdbReadHgncMapping <- function(path) {
    map <- read.delim(path, skip = 1L, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE)
    map <- map[grepl("ObjectDisplayForward", map[["GtP URL"]], fixed = TRUE), ]
    data.frame(targetId = as.character(map[["IUPHAR ID"]]),
              target_gene = map[["HGNC Symbol"]],
              stringsAsFactors = FALSE)
}

#' Read the interactions data from either route into one common shape
#'
#' \code{downloadGtoPdb()} can supply either the bulk interactions flat
#' file or a \code{/services/interactions} JSON response. Both are
#' mapped onto the same columns here so everything downstream - and the
#' \code{gtp_interactions} table itself - is identical either way. The
#' one column the REST route cannot supply is \code{pubmedId}: it
#' carries GtoPdb's internal reference IDs instead of the reference's
#' PubMed ID.
#' @keywords internal
#' @noRd
.gtoPdbReadInteractions <- function(path) {
    first <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
    if (grepl("^\\s*[[{]", first)) .gtoPdbFromJson(path)
    else .gtoPdbFromTsv(path)
}

#' @keywords internal
#' @noRd
.gtoPdbFromTsv <- function(path) {
    raw <- read.delim(path, skip = 1L, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      colClasses = "character")
    data.frame(
        targetId          = raw[["Target ID"]],
        targetName        = raw[["Target"]],
        targetSpecies     = .gtoPdbBlankToNa(raw[["Target Species"]]),
        primaryTarget     = raw[["Primary Target"]] == "true",
        ligandId          = raw[["Ligand ID"]],
        ligandName        = raw[["Ligand"]],
        type              = .gtoPdbBlankToNa(raw[["Type"]]),
        action            = .gtoPdbBlankToNa(raw[["Action"]]),
        affinity          = .gtoPdbAffinity(raw[["Affinity Low"]],
                                            raw[["Affinity Median"]],
                                            raw[["Affinity High"]]),
        affinityParameter = .gtoPdbBlankToNa(sub("^-$", "",
                                                 raw[["Affinity Units"]])),
        selectivity       = .gtoPdbBlankToNa(raw[["Selectivity"]]),
        pubmedId          = .gtoPdbBlankToNa(raw[["PubMed ID"]]),
        stringsAsFactors  = FALSE)
}

#' @keywords internal
#' @noRd
.gtoPdbFromJson <- function(path) {
    ix <- jsonlite::fromJSON(path)
    data.frame(
        targetId          = as.character(ix$targetId),
        targetName        = ix$targetName,
        targetSpecies     = ix$targetSpecies,
        primaryTarget     = as.logical(ix$primaryTarget),
        ligandId          = as.character(ix$ligandId),
        ligandName        = ix$ligandName,
        type              = ix$type,
        action            = ix$action,
        affinity          = ix$affinity,
        affinityParameter = ix$affinityParameter,
        selectivity       = ix$selectivity,
        pubmedId          = NA_character_,
        stringsAsFactors  = FALSE)
}


## ---------------------------------------------------------------------
## Build (or reuse) a local SQLite, versioned from the bulk files' own
## embedded release version - mirrors buildTtdDb()'s/
## buildBroadRepurposingHubDb()'s pattern of caching one ready-to-query
## SQLite file via BiocFileCache.
## ---------------------------------------------------------------------

#' Build (or fetch a cached) local SQLite database of GtoPdb drug-target
#' interactions
#'
#' Downloads GtoPdb's bulk interaction data and HGNC mapping (see
#' \code{\link{downloadGtoPdb}}), filters interactions to
#' \code{targetSpecies == "Human"}, left-joins in the HGNC gene symbol
#' per \code{targetId} (a target ID may map to more than one gene for
#' multi-subunit complexes - kept as separate rows, not collapsed), and
#' writes a single denormalized, indexed \code{gtp_interactions} table
#' to a local SQLite file - the same one-table shape as
#' \code{\link{buildTtdDb}}/\code{\link{buildBroadRepurposingHubDb}}. The
#' file is cached via BiocFileCache under a name that includes GtoPdb's
#' own release version (e.g. \code{gtopdb_2026.3.1.db}), so rebuilding is
#' a no-op until GtoPdb actually republishes. The package itself never
#' ships the full table - it is built into the caller's own local cache
#' the first time this is run.
#'
#' GtoPdb's bulk download files need no API key. To build from its
#' authenticated REST endpoint instead, pass \code{useServices = TRUE}
#' and an \code{apiKey}; see \code{\link{downloadGtoPdb}}.
#'
#' To build from files already on disk rather than downloading - the
#' small extract this package ships for its examples and vignette, or a
#' previously archived copy - pass \code{paths}.
#'
#' @param rerun logical(1); passed to \code{\link{downloadGtoPdb}}, and
#'   also controls whether an existing cached SQLite for the current
#'   version is reused (\code{FALSE}, the default here) or rebuilt
#'   (\code{TRUE}) - same rationale as \code{\link{buildTtdDb}}/
#'   \code{\link{buildBroadRepurposingHubDb}}: routine
#'   examples/tests/vignette chunks call this bare, and GtoPdb releases
#'   infrequently enough that defaulting to \code{TRUE} would just
#'   accumulate redundant cached copies.
#' @param config list as returned by \code{genConfig()}.
#' @param useServices logical(1); passed to
#'   \code{\link{downloadGtoPdb}}.
#' @param apiKey character(1); passed to \code{\link{downloadGtoPdb}}.
#' @param paths named list of local file paths with components
#'   \code{interactions} and \code{hgncMapping}, in the shape
#'   \code{\link{downloadGtoPdb}} returns. When supplied, nothing is
#'   downloaded and the database is built from these files into a
#'   temporary file whose path is returned, leaving the BiocFileCache
#'   copy of the full table untouched.
#' @return character(1) local file path to the SQLite database.
#' @examples
#' ## A small extract of the real GtoPdb interactions and HGNC mapping
#' ## files is shipped with the package, so the whole build-then-query
#' ## workflow runs here without downloading anything:
#' ex <- list(
#'     interactions = system.file("extdata",
#'         "gtopdb_interactions_example.tsv",
#'         package = "drugTargetInteractions"),
#'     hgncMapping = system.file("extdata",
#'         "gtopdb_hgnc_mapping_example.tsv",
#'         package = "drugTargetInteractions"))
#' dbPath <- buildGtoPdbDb(paths = ex)
#' gtoPdbTargetAnnot(
#'     list(molType = "protein", idType = "symbol", ids = "FGFR1"),
#'     dbPath)[1:3, c("ligandName", "type", "action", "affinity")]
#'
#' \donttest{
#'   ## The full table, downloaded and cached on first use:
#'   dbPath <- buildGtoPdbDb()
#'   dbPath
#' }
#' @seealso \code{\link{downloadGtoPdb}}, \code{\link{gtoPdbTargetAnnot}}
#' @export
buildGtoPdbDb <- function(rerun = FALSE, config = genConfig(),
                          useServices = FALSE,
                          apiKey = Sys.getenv("GTOPDB_API_KEY"),
                          paths = NULL) {
    supplied <- !is.null(paths)
    if (supplied) {
        if (!is.list(paths) ||
            !all(c("interactions", "hgncMapping") %in% names(paths)))
            stop("'paths' must be a named list with components ",
                 "'interactions' and 'hgncMapping', as returned by ",
                 "downloadGtoPdb().", call. = FALSE)
        want <- unlist(paths[c("interactions", "hgncMapping")])
        if (any(!file.exists(want)))
            stop("file(s) named in 'paths' do not exist: ",
                 paste(want[!file.exists(want)], collapse = ", "),
                 call. = FALSE)
    } else {
        paths <- downloadGtoPdb(rerun = rerun, config = config,
                                useServices = useServices, apiKey = apiKey)
    }

    ## A supplied-paths build is not the full table, so it must never be
    ## cached under - and so shadow - the real gtopdb_<version>.db.
    if (!supplied) {
        version <- .gtoPdbVersion(paths$hgncMapping)
        if (is.na(version)) version <- format(Sys.Date(), "%Y%m%d")
        dbName <- paste0("gtopdb_", version, ".db")
        if (!rerun) {
            existing <- tryCatch(.getCacheFile(dbName),
                                 error = function(e) NA_character_)
            if (length(existing) > 0 && !is.na(existing)) return(existing)
        }
    }

    ix <- .gtoPdbReadInteractions(paths$interactions)
    ix <- ix[!is.na(ix$targetSpecies) & ix$targetSpecies == "Human", ]
    ## Ligands whose target is not a GtoPdb object record carry no target
    ## identity at all, so they cannot form a (target, drug) edge - see
    ## this file's header.
    ix <- ix[nzchar(ix$targetId), ]
    ix$targetName <- .gtoPdbStripHtml(ix$targetName)
    ix$ligandName <- .gtoPdbStripHtml(ix$ligandName)

    hgncMap <- .gtoPdbReadHgncMapping(paths$hgncMapping)
    interactions <- merge(ix, hgncMap, by = "targetId", all.x = TRUE)
    interactions <- interactions[, c("target_gene", "targetId", "targetName",
                                     "targetSpecies", "primaryTarget", "ligandId",
                                     "ligandName", "type", "action", "affinity",
                                     "affinityParameter", "selectivity",
                                     "pubmedId")]
    colnames(interactions)[colnames(interactions) == "targetSpecies"] <- "species"

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "gtp_interactions", interactions, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_gtp_target ON gtp_interactions (target_gene)")
    dbExecute(con, "CREATE INDEX idx_gtp_ligand ON gtp_interactions (ligandName)")
    dbExecute(con, "CREATE INDEX idx_gtp_tid    ON gtp_interactions (targetId)")
    dbDisconnect(con)

    if (supplied) return(tmpDb)

    bfc <- .getCache()
    rid <- names(bfcadd(bfc, dbName, tmpDb, action = "copy"))
    file.remove(tmpDb)
    bfcrpath(bfc, rids = rid)
}

## ---------------------------------------------------------------------
## Bidirectional query, mirroring ttdTargetAnnot()'s/
## broadRepurposingHubAnnot()'s queryBy convention
## ---------------------------------------------------------------------

#' Query GtoPdb drug-target interactions bidirectionally
#'
#' Queries the local GtoPdb SQLite built by \code{\link{buildGtoPdbDb}},
#' using the same \code{queryBy = list(molType, idType, ids)} convention
#' as \code{\link{ttdTargetAnnot}}/\code{\link{broadRepurposingHubAnnot}}
#' - including its \code{QueryIDs} column: every row of the result is
#' tagged with the original query token it matched, and query IDs that
#' returned no rows still appear as a single row with all other fields
#' \code{NA}.
#' \itemize{
#'   \item \code{molType = "protein"} (or \code{"gene"}), \code{idType}
#'     one of \code{"symbol"} (HGNC gene symbol) or
#'     \code{"gtp_target_id"} (GtoPdb's own numeric target ID) -> target
#'     -> drug.
#'   \item \code{molType = "cmp"}, \code{idType} one of \code{"name"}
#'     (ligand name) or \code{"gtp_ligand_id"} (GtoPdb's own numeric
#'     ligand ID) -> drug -> target.
#' }
#' \code{"symbol"}/\code{"name"} lookups are case-insensitive; the
#' \code{"gtp_target_id"}/\code{"gtp_ligand_id"} ID lookups are
#' exact-match.
#'
#' @param queryBy list with components \code{molType}, \code{idType},
#'   \code{ids} (character vector).
#' @param gtoPdbDbPath character(1) path to the GtoPdb SQLite, e.g. from
#'   \code{\link{buildGtoPdbDb}}.
#' @param fields \code{"core"} (default) or \code{"all"} - both return
#'   every \code{gtp_interactions} column, since (unlike the REST-backed
#'   sources) there is no larger raw payload to opt into - or a character
#'   vector of column names to keep (\code{QueryIDs} is always retained).
#'   See \code{\link{listDrugTargetFields}}.
#' @return A \code{data.frame} with columns \code{QueryIDs},
#'   \code{target_gene}, \code{targetId}, \code{targetName}, \code{species},
#'   \code{primaryTarget}, \code{ligandId}, \code{ligandName}, \code{type},
#'   \code{action}, \code{affinity}, \code{affinityParameter},
#'   \code{selectivity}, \code{pubmedId} (with the default
#'   \code{fields = "core"}), or a subset when \code{fields} requests
#'   specific columns.
#' @examples
#' ## A small extract of GtoPdb's two bulk files ships with the package,
#' ## so the whole build-then-query path runs without downloading.
#' ex <- list(
#'     interactions = system.file("extdata", "gtopdb_interactions_example.tsv",
#'                                package = "drugTargetInteractions"),
#'     hgncMapping = system.file("extdata", "gtopdb_hgnc_mapping_example.tsv",
#'                               package = "drugTargetInteractions"))
#' dbPath <- buildGtoPdbDb(paths = ex)
#'
#' ## Target -> drug
#' head(gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol",
#'                             ids = "FGFR1"), dbPath)[, c("ligandName", "action")], 3)
#'
#' ## Drug -> target
#' gtoPdbTargetAnnot(list(molType = "cmp", idType = "name",
#'                        ids = "pemigatinib"), dbPath)[, c("target_gene", "affinity")]
#'
#' \donttest{
#'   dbPath <- buildGtoPdbDb()
#'   gtoPdbTargetAnnot(list(molType = "protein", idType = "symbol",
#'                          ids = c("FGFR1", "IL1B")), dbPath)
#'   gtoPdbTargetAnnot(list(molType = "cmp", idType = "name",
#'                          ids = "pemigatinib"), dbPath)
#' }
#' @seealso \code{\link{buildGtoPdbDb}}, \code{\link{ttdTargetAnnot}},
#'   \code{\link{broadRepurposingHubAnnot}}, \code{\link{listDrugTargetFields}}
#' @export
gtoPdbTargetAnnot <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                              gtoPdbDbPath, fields = "core") {
    if (any(names(queryBy) != c("molType", "idType", "ids"))) {
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

    col <- if (queryBy$molType %in% c("protein", "gene")) {
        switch(queryBy$idType,
               symbol        = "target_gene",
               gtp_target_id = "targetId",
               stop("idType for molType='protein'/'gene' must be one of: ",
                    "'symbol', 'gtp_target_id'"))
    } else if (queryBy$molType == "cmp") {
        switch(queryBy$idType,
               name          = "ligandName",
               gtp_ligand_id = "ligandId",
               stop("idType for molType='cmp' must be one of: ",
                    "'name', 'gtp_ligand_id'"))
    } else {
        stop("molType must be 'protein'/'gene' or 'cmp'")
    }

    caseInsensitive <- queryBy$idType %in% c("symbol", "name")
    ids <- if (caseInsensitive) toupper(queryBy$ids) else queryBy$ids
    idvec <- paste0("('", paste(gsub("'", "''", ids, fixed = TRUE), collapse = "', '"), "')")
    colExpr <- if (caseInsensitive) paste0("UPPER(", col, ")") else col

    con <- dbConnect(SQLite(), gtoPdbDbPath)
    on.exit(dbDisconnect(con))
    query <- paste0("SELECT * FROM gtp_interactions WHERE ", colExpr, " IN ", idvec)
    resultDF <- dbGetQuery(con, query)

    ## Tag every row with the original query token it matched, and NA-pad
    ## any query ID that matched nothing - mirrors ttdTargetAnnot()'s/
    ## broadRepurposingHubAnnot()'s QueryIDs convention exactly (same
    ## Inf-index trick: an unmatched ID gets rowid Inf, and indexing a
    ## data.frame with Inf yields a row of NAs).
    cmpFun <- if (caseInsensitive) toupper else identity
    index_list <- lapply(
        queryBy$ids,
        function(x) which(cmpFun(resultDF[[col]]) %in% cmpFun(x))
    )
    names(index_list) <- queryBy$ids
    index_list[vapply(index_list, length, integer(1)) == 0] <- Inf
    index_df <- data.frame(
        ids = rep(names(index_list), vapply(index_list, length, integer(1))),
        rowids = unlist(index_list)
    )
    out <- data.frame(
        QueryIDs = index_df[, 1],
        resultDF[as.numeric(index_df$rowids), ],
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    .dtiSelectFields(out, fields, .dtiGtoPdbAllCols)
}

#' Documented column list for \code{listDrugTargetFields("gtopdb")}
#'
#' Unlike the REST-backed sources' \code{fields = "all"} (ChEMBL,
#' PubChem, DGIdb, Open Targets), \code{gtp_interactions} is a single
#' flat local SQLite table built entirely by \code{\link{buildGtoPdbDb}}
#' (see there), so this is an exact list, not a best-effort one, and
#' \code{fields = "core"} and \code{fields = "all"} are equivalent for
#' \code{\link{gtoPdbTargetAnnot}}.
#' @keywords internal
#' @noRd
.dtiGtoPdbAllCols <- c("QueryIDs", "target_gene", "targetId", "targetName",
                       "species", "primaryTarget", "ligandId", "ligandName",
                       "type", "action", "affinity", "affinityParameter",
                       "selectivity", "pubmedId")
