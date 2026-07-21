## =====================================================================
##  gtoPdbAccess.R
##  IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) integration for the
##  drugTargetInteractions Bioconductor package.
##
##  Unlike TTD/the Broad Repurposing Hub (no API at all), GtoPdb offers a
##  documented REST API (https://www.guidetopharmacology.org/webServices.jsp)
##  - but its per-target/per-ligand endpoints are one-resource-per-call,
##  with no bulk endpoint for target ID -> gene symbol mapping (only a
##  per-target databaseLinks call, which would mean 3000+ requests for a
##  genome-wide build). This module instead uses GtoPdb's *bulk* REST
##  endpoint for the actual interaction data (a single call returns the
##  whole ~24,000-row interactions table as JSON - live-confirmed
##  2026-07-19, no pagination), joined locally against the small bulk
##  flat file GtoPdb separately publishes for target ID <-> HGNC symbol
##  mapping (GtP_to_HGNC_mapping.tsv, ~3,800 rows). Both are downloaded
##  once (cached via the package's existing BiocFileCache helpers) and
##  built into a local SQLite - the same build-once-query-many pattern as
##  ttdAccess.R/broadRepurposingHubAccess.R, just sourced from REST + one
##  small file instead of pure flat files.
##
##  LICENSE / distribution posture: GtoPdb data is licensed under the
##  Open Data Commons Open Database License (ODbL) with contents under
##  CC-BY-SA 4.0 - clear terms that explicitly permit redistribution
##  (with attribution/share-alike), unlike TTD's ambiguous or the Broad
##  Repurposing Hub's non-commercial-only wording. This module still
##  defaults to the same code-only, build-locally-at-call-time posture as
##  the other local-SQLite sources for consistency - the package itself
##  ships no GtoPdb data - but redistribution would not be a licensing
##  problem here if that default is ever revisited.
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
##  Interactions are filtered to targetSpecies == "Human" - the rest of
##  this package is implicitly human-gene-centric (HGNC anchor, human
##  UniProt accessions), and GtP_to_HGNC_mapping.tsv itself mixes
##  human/mouse/rat gene symbols in its "HGNC Symbol" column for
##  cross-species target entries, so species-filtering the interactions
##  side first avoids any risk of a cross-species symbol collision.
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for GtoPdb's bulk REST endpoint and HGNC mapping file
#' @keywords internal
.gtoPdbEndpoints <- function() {
    list(
        interactions = "https://www.guidetopharmacology.org/services/interactions",
        hgncMapping  = "https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.tsv"
    )
}

#' Strip HTML markup (GtoPdb target/ligand names embed \\code{<sub>}/
#' \\code{<sup>}/\\code{<i>}/\\code{<small>} tags, e.g.
#' \\code{"5-HT<sub>1A</sub> receptor"}) down to plain text.
#' @keywords internal
.gtoPdbStripHtml <- function(x) {
    gsub("<[^>]+>", "", x)
}

#' Extract GtoPdb's own release version from GtP_to_HGNC_mapping.tsv's
#' "# GtoPdb Version: X.Y - published: YYYY-MM-DD" first-line comment.
#' @keywords internal
.gtoPdbVersion <- function(path) {
    hdr <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
    hit <- regmatches(hdr, regexpr("Version:\\s*[0-9.]+", hdr))
    if (length(hit) == 0L) return(NA_character_)
    trimws(sub("Version:\\s*", "", hit))
}


## ---------------------------------------------------------------------
## Download raw data (cached via the package's existing BiocFileCache
## helpers)
## ---------------------------------------------------------------------

#' Download GtoPdb's bulk interactions (REST) and HGNC mapping (flat
#' file) data
#'
#' Downloads (or reuses previously cached copies of) GtoPdb's bulk
#' \code{/services/interactions} REST response and its
#' \code{GtP_to_HGNC_mapping.tsv} flat file via the package's existing
#' \code{.downloadFile()}/BiocFileCache infrastructure - the same
#' mechanism \code{downloadTTD()}/\code{downloadBroadRepurposingHub()}
#' use. Files land in the local BiocFileCache (see \code{.getCache()},
#' \code{rappdirs::user_cache_dir(appname = "drugTargetInteractions")}) -
#' no GtoPdb data is bundled with or downloaded by the package itself
#' until this is called explicitly.
#'
#' @param rerun logical(1); if \code{TRUE} (default), check for updates
#'   and (re)download as needed; if \code{FALSE}, use whatever is
#'   already in the local cache without checking upstream.
#' @param config list as returned by \code{genConfig()}; unused here but
#'   accepted for consistency with the rest of the package's API.
#' @return A named list of local file paths: \code{interactions},
#'   \code{hgncMapping}.
#' @examples
#' \donttest{
#'   paths <- downloadGtoPdb()
#'   paths
#' }
#' @seealso \code{\link{buildGtoPdbDb}}
#' @export
downloadGtoPdb <- function(rerun = TRUE, config = genConfig()) {
    ep <- .gtoPdbEndpoints()
    fnames <- c(interactions = "gtp_interactions.json",
                hgncMapping  = "GtP_to_HGNC_mapping.tsv")
    paths <- vector("list", length(ep))
    names(paths) <- names(ep)
    for (nm in names(ep)) {
        if (rerun) {
            paths[[nm]] <- .downloadFile(ep[[nm]], fnames[[nm]])
        } else {
            paths[[nm]] <- .getCacheFile(fnames[[nm]])
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
.gtoPdbReadHgncMapping <- function(path) {
    map <- read.delim(path, skip = 1L, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE, check.names = FALSE)
    map <- map[grepl("ObjectDisplayForward", map[["GtP URL"]], fixed = TRUE), ]
    data.frame(targetId = map[["IUPHAR ID"]], target_gene = map[["HGNC Symbol"]],
              stringsAsFactors = FALSE)
}


## ---------------------------------------------------------------------
## Build (or reuse) a local SQLite, versioned from GtP_to_HGNC_mapping's
## own embedded release version - mirrors buildTtdDb()'s/
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
#' own release version (e.g. \code{gtopdb_2026.2.db}), so rebuilding is a
#' no-op until GtoPdb actually republishes. As with the other local-SQLite
#' sources, the package itself never ships or redistributes this file -
#' it is built into the caller's own local cache the first time this is
#' run (see the file header for why that is a conservative default here,
#' not a licensing requirement).
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
#' @return character(1) local file path to the SQLite database.
#' @examples
#' \donttest{
#'   dbPath <- buildGtoPdbDb()
#'   dbPath
#' }
#' @seealso \code{\link{downloadGtoPdb}}, \code{\link{gtoPdbTargetAnnot}}
#' @export
buildGtoPdbDb <- function(rerun = FALSE, config = genConfig()) {
    paths <- downloadGtoPdb(rerun = rerun, config = config)

    version <- .gtoPdbVersion(paths$hgncMapping)
    if (is.na(version)) version <- format(Sys.Date(), "%Y%m%d")
    dbName <- paste0("gtopdb_", version, ".db")

    if (!rerun) {
        existing <- tryCatch(.getCacheFile(dbName), error = function(e) NA_character_)
        if (length(existing) > 0 && !is.na(existing)) return(existing)
    }

    ix <- jsonlite::fromJSON(paths$interactions)
    ix <- ix[!is.na(ix$targetSpecies) & ix$targetSpecies == "Human", ]
    ix$targetName <- .gtoPdbStripHtml(ix$targetName)
    ix$ligandName <- .gtoPdbStripHtml(ix$ligandName)
    ## refIds/dataPointInteractionIds come back from jsonlite as
    ## list-columns (some rows have multiple values) - SQLite can only
    ## store atomic columns, so flatten refIds to a single comma-joined
    ## string; dataPointInteractionIds is dropped entirely (internal
    ## bookkeeping ID, not useful downstream).
    ix$refIds <- vapply(ix$refIds, paste, character(1), collapse = ", ")

    hgncMap <- .gtoPdbReadHgncMapping(paths$hgncMapping)
    interactions <- merge(ix, hgncMap, by = "targetId", all.x = TRUE)
    interactions <- interactions[, c("target_gene", "targetId", "targetName",
                                     "targetSpecies", "primaryTarget", "ligandId",
                                     "ligandName", "type", "action", "affinity",
                                     "affinityParameter", "selectivity", "refIds")]
    colnames(interactions)[colnames(interactions) == "targetSpecies"] <- "species"

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "gtp_interactions", interactions, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_gtp_target ON gtp_interactions (target_gene)")
    dbExecute(con, "CREATE INDEX idx_gtp_ligand ON gtp_interactions (ligandName)")
    dbExecute(con, "CREATE INDEX idx_gtp_tid    ON gtp_interactions (targetId)")
    dbDisconnect(con)

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
#'   \code{selectivity}, \code{refIds} (with the default
#'   \code{fields = "core"}), or a subset when \code{fields} requests
#'   specific columns.
#' @examples
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
.dtiGtoPdbAllCols <- c("QueryIDs", "target_gene", "targetId", "targetName",
                       "species", "primaryTarget", "ligandId", "ligandName",
                       "type", "action", "affinity", "affinityParameter",
                       "selectivity", "refIds")
