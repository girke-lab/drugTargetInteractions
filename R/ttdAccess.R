## =====================================================================
##  ttdAccess.R
##  TTD (Therapeutic Targets Database) integration for the
##  drugTargetInteractions Bioconductor package.
##
##  TTD has NO programmatic API; bulk flat files are served from
##  https://ttd.idrblab.cn/files/download/<FILENAME>. This module
##  downloads the needed files once (cached via the package's existing
##  BiocFileCache helpers - .getCache()/.getCacheFile()/.downloadFile(),
##  defined in drugTargetAnnotations_Fct.R), parses them into a
##  normalized table, and writes it to a local SQLite database. This
##  mirrors downloadChemblDb()'s local-SQL pattern rather than
##  apiAccess.R's live-API pattern, since TTD data has to be downloaded
##  and parsed before it is queryable - kept in its own file for the same
##  reason downloadChemblDb()/drugTargetAnnot() and this module are
##  conceptually paired: build/cache a local database once, then query it
##  repeatedly, as opposed to a live per-call REST/GraphQL round trip.
##
##  LICENSE / distribution posture: TTD states data is "freely available
##  ... for academic use" but publishes no explicit redistribution
##  license (no CC-BY / CC0). This module only ever downloads TTD's own
##  files into the *caller's* local BiocFileCache at call time and builds
##  a local SQLite there - the package itself ships no TTD data file and
##  never redistributes TTD content, matching the posture already used
##  for ChEMBL's downloadChemblDb(). Until/unless TTD's authors grant
##  redistribution permission, this stays fetch-and-build-locally only.
##
##  Files used (see .ttdEndpoints()):
##    P1-01-TTD_target_download.txt   target records incl. per-target
##                                     DRUGINFO lines (TargetID, DrugID,
##                                     DrugName, Highest Clinical Status)
##    P1-02-TTD_drug_download.txt     drug records incl. DRUGSMIL (SMILES)
##    P1-07-Drug-TargetMapping.xlsx   TargetID x DrugID x Highest_status x
##                                     MOA - the actual interaction edges
##
##  Every TTD flat file embeds its own release version on line 3 of its
##  header, e.g. "Version 10.1.01 (2024.01.10)" - this is TTD's real
##  upstream version (their equivalent of a ChEMBL release number) and is
##  used to name/cache the derived SQLite file, NOT an arbitrary download
##  timestamp.
##
##  Suggested DESCRIPTION addition: Imports: readxl
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for TTD flat files
#' @keywords internal
.ttdEndpoints <- function() {
    list(
        base    = "https://ttd.idrblab.cn/files/download",
        targets = "P1-01-TTD_target_download.txt",
        drugs   = "P1-02-TTD_drug_download.txt",
        mapping = "P1-07-Drug-TargetMapping.xlsx"
    )
}

#' Extract TTD's own release version from a downloaded flat file header
#'
#' TTD stamps every flat file with "Version X.Y.Z (YYYY.MM.DD)" on line 3
#' of its header block. Used to name the derived SQLite so re-running the
#' build is a no-op until TTD actually bumps this string.
#'
#' @param path character(1) path to a downloaded TTD flat file.
#' @return character(1) version string (e.g. \code{"10.1.01"}), or
#'   \code{NA} if not found.
#' @keywords internal
.ttdVersion <- function(path) {
    hdr <- readLines(path, n = 5L, warn = FALSE, encoding = "UTF-8")
    hit <- grep("^Version ", hdr, value = TRUE)
    if (length(hit) == 0L) return(NA_character_)
    sub("^Version ([0-9.]+).*$", "\\1", hit[1])
}


## ---------------------------------------------------------------------
## Download raw TTD flat files (cached via the package's existing
## BiocFileCache helpers)
## ---------------------------------------------------------------------

#' Download the TTD flat files needed for target-drug mapping
#'
#' Downloads (or reuses previously cached copies of) the TTD target,
#' drug and drug-target-mapping flat files via the package's existing
#' \code{.downloadFile()}/BiocFileCache infrastructure - the same
#' mechanism \code{downloadChemblDb()} and \code{downloadUniChem()} use.
#' Files land in the local BiocFileCache (see \code{.getCache()},
#' \code{rappdirs::user_cache_dir(appname = "drugTargetInteractions")|}) -
#' no TTD data is bundled with or downloaded by the package itself until
#' this is called explicitly.
#'
#' @param rerun logical(1); if \code{TRUE} (default), check for updates
#'   and (re)download as needed; if \code{FALSE}, use whatever is
#'   already in the local cache without checking upstream.
#' @param config list as returned by \code{genConfig()}; unused here but
#'   accepted for consistency with the rest of the package's API.
#' @return A named list of local file paths: \code{targets}, \code{drugs},
#'   \code{mapping}.
#' @examples
#' \donttest{
#'   paths <- downloadTTD()
#'   paths
#' }
#' @seealso \code{\link{buildTtdDb}}
#' @export
downloadTTD <- function(rerun = TRUE, config = genConfig()) {
    ep    <- .ttdEndpoints()
    files <- ep[c("targets", "drugs", "mapping")]
    paths <- vector("list", length(files))
    names(paths) <- names(files)
    for (nm in names(files)) {
        fname <- files[[nm]]
        if (rerun) {
            paths[[nm]] <- .downloadFile(paste0(ep$base, "/", fname), fname)
        } else {
            paths[[nm]] <- .getCacheFile(fname)
        }
    }
    paths
}


## ---------------------------------------------------------------------
## Parse the raw TTD flat files
## ---------------------------------------------------------------------

## P1-01/P1-02 are line-oriented "ID <TAB> FIELD <TAB> VALUE [...]"
## dumps, not plain tables (one record spans many lines); IDs always
## start with "T" (targets) or "D" (drugs), which is used to strip the
## header/legend block without needing a stateful line scan.

#' Long (id, field, value) table for one TTD flat file
#' @keywords internal
.ttdReadLong <- function(path, idPrefix) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    parts <- strsplit(lines, "\t", fixed = TRUE)
    len   <- lengths(parts)
    keep  <- len >= 3L & startsWith(lines, idPrefix)
    parts <- parts[keep]
    data.frame(
        id    = vapply(parts, `[`, character(1), 1L),
        field = vapply(parts, `[`, character(1), 2L),
        value = vapply(parts, `[`, character(1), 3L),
        stringsAsFactors = FALSE
    )
}

#' First-value-wins lookup for a single field of a long table
#' @keywords internal
.ttdFieldLookup <- function(long, field) {
    sub <- long[long$field == field, c("id", "value")]
    sub[!duplicated(sub$id), ]
}

#' Parse P1-01 into one row per TargetID (GeneName, Uniprot, TargetType)
#' @keywords internal
.ttdParseTargets <- function(path) {
    long <- .ttdReadLong(path, "T")
    gene <- .ttdFieldLookup(long, "GENENAME")
    up   <- .ttdFieldLookup(long, "UNIPROID")
    typ  <- .ttdFieldLookup(long, "TARGTYPE")
    out  <- data.frame(TargetID = unique(long$id), stringsAsFactors = FALSE)
    out$GeneName   <- gene$value[match(out$TargetID, gene$id)]
    out$Uniprot    <- up$value[match(out$TargetID, up$id)]
    out$TargetType <- typ$value[match(out$TargetID, typ$id)]
    out
}

#' Drug-name lookup from P1-01's per-target DRUGINFO lines
#'
#' DRUGINFO lines carry (TargetID, "DRUGINFO", DrugID, DrugName, Status);
#' this collapses them to a global DrugID -> DrugName lookup (first name
#' wins), independent of which target(s) the drug maps to.
#' @keywords internal
.ttdParseTargetDrugNames <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- grep("\tDRUGINFO\t", lines, value = TRUE, fixed = TRUE)
    parts <- strsplit(lines, "\t", fixed = TRUE)
    parts <- parts[lengths(parts) >= 5L]
    out <- data.frame(
        DrugID   = vapply(parts, `[`, character(1), 3L),
        DrugName = vapply(parts, `[`, character(1), 4L),
        stringsAsFactors = FALSE
    )
    out[!duplicated(out$DrugID), ]
}

#' Parse P1-02's DRUGSMIL (canonical SMILES) field into a lookup
#' @keywords internal
.ttdParseDrugSmiles <- function(path) {
    long <- .ttdReadLong(path, "D")
    sm <- .ttdFieldLookup(long, "DRUGSMIL")
    colnames(sm) <- c("DrugID", "Smiles")
    sm
}


## ---------------------------------------------------------------------
## Build (or reuse) a local TTD SQLite, versioned from the raw files'
## own embedded release version - mirrors downloadChemblDb()'s pattern
## of caching one ready-to-query SQLite file via BiocFileCache.
## ---------------------------------------------------------------------

#' Build (or fetch a cached) local SQLite database of TTD interactions
#'
#' Downloads the TTD flat files (see \code{\link{downloadTTD}}), parses
#' them, and writes a single denormalized \code{ttd_interactions} table
#' (TTD's target-drug relationship is a simple 1-hop mapping, unlike
#' ChEMBL's multi-table mechanism/activity schema, so one flat, indexed
#' table is the natural queryable shape) to a local SQLite file. The
#' file is cached via BiocFileCache under a name that includes TTD's own
#' release version (e.g. \code{ttd_10.1.01.db}), so rebuilding is a
#' no-op until TTD actually publishes a new version. As with
#' \code{\link{downloadChemblDb}}, the package itself never ships or
#' redistributes this file - it is built into the caller's own local
#' cache the first time this is run.
#'
#' @param rerun logical(1); passed to \code{\link{downloadTTD}}, and
#'   also controls whether an existing cached SQLite for the current
#'   version is reused (\code{FALSE}, the default here) or rebuilt
#'   (\code{TRUE}). Deliberately defaults to \code{FALSE}, unlike most
#'   other \code{download*()}/\code{build*()} functions in this package
#'   (which default to \code{rerun = TRUE}): TTD publishes a new release
#'   only rarely, and a bare \code{buildTtdDb()} call is routine in
#'   examples, tests, and vignette chunks - defaulting to \code{TRUE}
#'   there meant every such run silently rebuilt and re-cached the
#'   database from scratch, accumulating redundant copies over time
#'   with no benefit. Pass \code{rerun = TRUE} explicitly to force
#'   a fresh check for a new TTD release.
#' @param config list as returned by \code{genConfig()}.
#' @return character(1) local file path to the SQLite database.
#' @examples
#' \donttest{
#'   dbPath <- buildTtdDb()
#'   dbPath
#' }
#' @seealso \code{\link{downloadTTD}}, \code{\link{ttdTargetAnnot}}
#' @export
buildTtdDb <- function(rerun = FALSE, config = genConfig()) {
    paths   <- downloadTTD(rerun = rerun, config = config)
    version <- .ttdVersion(paths$targets)
    if (is.na(version)) version <- format(Sys.Date(), "%Y%m%d")
    dbName  <- paste0("ttd_", version, ".db")

    if (!rerun) {
        existing <- tryCatch(.getCacheFile(dbName), error = function(e) NA_character_)
        if (length(existing) > 0 && !is.na(existing)) return(existing)
    }

    targets     <- .ttdParseTargets(paths$targets)
    drugNames   <- .ttdParseTargetDrugNames(paths$targets)
    drugSmiles  <- .ttdParseDrugSmiles(paths$drugs)
    mapping     <- as.data.frame(readxl::read_excel(paths$mapping))

    interactions <- merge(mapping, targets, by = "TargetID", all.x = TRUE)
    interactions <- merge(interactions, drugNames, by = "DrugID", all.x = TRUE)
    interactions <- merge(interactions, drugSmiles, by = "DrugID", all.x = TRUE)
    interactions <- interactions[, c("TargetID", "GeneName", "Uniprot", "TargetType",
                                     "DrugID", "DrugName", "Smiles",
                                     "Highest_status", "MOA")]

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "ttd_interactions", interactions, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_ttd_target   ON ttd_interactions (TargetID)")
    dbExecute(con, "CREATE INDEX idx_ttd_gene     ON ttd_interactions (GeneName)")
    dbExecute(con, "CREATE INDEX idx_ttd_drug     ON ttd_interactions (DrugID)")
    dbExecute(con, "CREATE INDEX idx_ttd_drugname ON ttd_interactions (DrugName)")
    dbDisconnect(con)

    bfc <- .getCache()
    rid <- names(bfcadd(bfc, dbName, tmpDb, action = "copy"))
    file.remove(tmpDb)
    bfcrpath(bfc, rids = rid)
}


## ---------------------------------------------------------------------
## Bidirectional query, mirroring drugTargetAnnot()'s queryBy convention
## ---------------------------------------------------------------------

#' Query TTD target-drug interactions bidirectionally
#'
#' Queries the local TTD SQLite built by \code{\link{buildTtdDb}}, using
#' the same \code{queryBy = list(molType, idType, ids)} convention as the
#' package's ChEMBL-SQL \code{\link{drugTargetAnnot}} - including its
#' \code{QueryIDs} column: every row of the result is tagged with the
#' original query token it matched, and query IDs that returned no rows
#' still appear as a single row with all other fields \code{NA}, so
#' callers can always confirm which of their input IDs were resolved
#' (unlike the standalone \code{R_Py_code/ttdAccess.R} reference this was
#' ported from, which returned only the matched rows with no such
#' bookkeeping).
#' \itemize{
#'   \item \code{molType = "protein"}, \code{idType} one of
#'     \code{"symbol"}, \code{"uniprot"}, \code{"ttd_target_id"} ->
#'     target -> drug
#'   \item \code{molType = "cmp"}, \code{idType} one of \code{"name"},
#'     \code{"ttd_drug_id"} -> drug -> target
#' }
#' Symbol/name lookups are case-insensitive; ID lookups
#' (\code{"uniprot"}, \code{"ttd_target_id"}, \code{"ttd_drug_id"}) are
#' exact-match. Note: TTD's \code{"uniprot"} field is the UniProt
#' *mnemonic entry name* (e.g. \code{"FGFR1_HUMAN"}), NOT the accession
#' number (e.g. \code{"P11362"}) used by the package's ChEMBL-SQL side
#' (\code{\link{drugTargetAnnot}}) - a real cross-source naming
#' inconsistency, not a bug; convert accession -> mnemonic upstream if
#' chaining the two.
#'
#' @param queryBy list with components \code{molType}, \code{idType},
#'   \code{ids} (character vector).
#' @param ttdDbPath character(1) path to the TTD SQLite, e.g. from
#'   \code{\link{buildTtdDb}}.
#' @return A \code{data.frame} with columns \code{QueryIDs},
#'   \code{TargetID}, \code{GeneName}, \code{Uniprot}, \code{TargetType},
#'   \code{DrugID}, \code{DrugName}, \code{Smiles}, \code{Highest_status},
#'   \code{MOA}.
#' @examples
#' \donttest{
#'   dbPath <- buildTtdDb()
#'   ttdTargetAnnot(list(molType = "protein", idType = "symbol",
#'                        ids = c("FGFR1", "IL1B")), dbPath)
#'   ttdTargetAnnot(list(molType = "cmp", idType = "name",
#'                        ids = "Pemigatinib"), dbPath)
#' }
#' @seealso \code{\link{buildTtdDb}}, \code{\link{drugTargetAnnot}}
#' @export
ttdTargetAnnot <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                           ttdDbPath) {
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

    col <- if (queryBy$molType == "protein") {
        switch(queryBy$idType,
               symbol        = "GeneName",
               uniprot       = "Uniprot",
               ttd_target_id = "TargetID",
               stop("idType for molType='protein' must be one of: ",
                    "'symbol', 'uniprot', 'ttd_target_id'"))
    } else if (queryBy$molType == "cmp") {
        switch(queryBy$idType,
               name        = "DrugName",
               ttd_drug_id = "DrugID",
               stop("idType for molType='cmp' must be one of: ",
                    "'name', 'ttd_drug_id'"))
    } else {
        stop("molType must be 'protein' or 'cmp'")
    }

    caseInsensitive <- queryBy$idType %in% c("symbol", "name")
    ids <- if (caseInsensitive) toupper(queryBy$ids) else queryBy$ids
    idvec <- paste0("('", paste(gsub("'", "''", ids, fixed = TRUE), collapse = "', '"), "')")
    colExpr <- if (caseInsensitive) paste0("UPPER(", col, ")") else col

    con <- dbConnect(SQLite(), ttdDbPath)
    on.exit(dbDisconnect(con))
    query <- paste0("SELECT * FROM ttd_interactions WHERE ", colExpr, " IN ", idvec)
    resultDF <- dbGetQuery(con, query)

    ## Tag every row with the original query token it matched, and NA-pad
    ## any query ID that matched nothing - mirrors drugTargetAnnot()'s
    ## QueryIDs convention exactly (same Inf-index trick: an unmatched ID
    ## gets rowid Inf, and indexing a data.frame with Inf yields a row of
    ## NAs), just joining on `col` (the idType's resolved column name)
    ## rather than `queryBy$idType` itself, since TTD's idType vocabulary
    ## ("symbol", "uniprot", ...) is descriptive, not a literal column name.
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
    out
}
