## =====================================================================
##  broadRepurposingHubAccess.R
##  Broad Institute Drug Repurposing Hub (CLUE) integration for the
##  drugTargetInteractions Bioconductor package.
##
##  The Repurposing Hub has no programmatic API; two flat TSV files are
##  served directly (drug-level annotation, keyed by `pert_iname`, and
##  sample-level annotation, keyed by the physical vial/batch ID
##  `broad_id`). This mirrors ttdAccess.R's pattern exactly: download the
##  flat files once (cached via the package's existing BiocFileCache
##  helpers - .getCache()/.getCacheFile()/.downloadFile(), defined in
##  drugTargetAnnotations_Fct.R), parse and join them, and write a single
##  denormalized, indexed table to a local SQLite - build once, query
##  many times.
##
##  LICENSE / distribution posture: both files explicitly state ("!Restrictions"
##  header line) "The Drug Repurposing Hub data is provided for
##  non-commercial use only" - an explicit restriction, stricter than
##  TTD's merely-ambiguous "freely available for academic use" wording.
##  As with downloadTTD()/buildTtdDb(), this module only ever downloads
##  the Hub's own files into the *caller's* local BiocFileCache at call
##  time and builds a local SQLite there - the package itself ships no
##  Repurposing Hub data and never redistributes it: code only, never
##  data.
##
##  Files used (see .brhEndpoints()) - both are TSVs with a 9-line
##  "!Key\tValue" metadata block before the real header row:
##    repo-drug-annotation-<date>.txt    one row per drug (pert_iname):
##                                        clinical_phase, moa, target
##                                        (" | "-separated gene symbols,
##                                        may be empty), disease_area,
##                                        indication
##    repo-sample-annotation-<date>.txt  one row per physical sample
##                                        (broad_id): links back to
##                                        pert_iname, plus vendor/purity/
##                                        structure (smiles, InChIKey,
##                                        pubchem_cid)
##
##  The date embedded in each URL/filename is NOT a reliable version
##  indicator - both files' own "!File_date" metadata line (the real
##  content version) was observed to read the same current date despite
##  the URLs themselves carrying stale 2020/2024 filename suffixes. As
##  with TTD's embedded "Version X.Y.Z" line, "!File_date" (not the URL)
##  is what's used to name/cache the derived SQLite. Because the URLs
##  themselves are otherwise fixed/hardcoded, a future Broad Hub release
##  that changes the filename suffix will need .brhEndpoints() updated by
##  hand - there is no discovery mechanism for the current filename.
##
##  Known operational caveat, unrelated to this code: as of this writing
##  repo-hub.broadinstitute.org's TLS certificate chain is served without
##  its InCommon intermediate certificate, which can cause certificate
##  verification failures on some machines/OSes whose trust store hasn't
##  independently cached that intermediate (most browsers chase the AIA
##  "CA Issuers" URL automatically; curl/libcurl/R's default download
##  methods generally do not). This is a server-side misconfiguration -
##  do not "fix" it by disabling certificate verification here.
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for the Broad Repurposing Hub flat files
#' @keywords internal
.brhEndpoints <- function() {
    list(
        drug   = "https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20200324.txt",
        sample = "https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20240610.txt"
    )
}

#' Parse a Repurposing Hub flat file: skip the "!Key\\tValue" metadata
#' block, read the real header + data, and stash the file's own
#' "!File_date" as an attribute (the real version indicator - see the
#' file header for why the URL's own embedded date is not reliable).
#'
#' Uses default (non-disabled) quote handling, like the package's HGNC
#' parser (\code{getHgncGeneTable()}, genomeWideAnnot.R) - these files
#' rely on ordinary CSV-style quoting for fields containing commas (e.g.
#' \code{moa} values), and disabling it would leave stray literal quote
#' characters in the parsed values.
#' @keywords internal
.brhReadTable <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    meta  <- grep("^!", lines)
    hdrIdx <- max(meta) + 1L
    df <- read.delim(text = paste(lines[hdrIdx:length(lines)], collapse = "\n"),
                      sep = "\t", quote = "\"", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = "")
    fdLine <- grep("^!File_date", lines, value = TRUE)
    attr(df, "fileDate") <- if (length(fdLine)) trimws(strsplit(fdLine[1], "\t")[[1]][2]) else NA_character_
    df
}


## ---------------------------------------------------------------------
## Download raw flat files (cached via the package's existing
## BiocFileCache helpers)
## ---------------------------------------------------------------------

#' Download the Broad Repurposing Hub drug/sample annotation flat files
#'
#' Downloads (or reuses previously cached copies of) the two Repurposing
#' Hub flat files via the package's existing \code{.downloadFile()}/
#' BiocFileCache infrastructure - the same mechanism \code{downloadTTD()}
#' uses. Files land in the local BiocFileCache (see \code{.getCache()},
#' \code{rappdirs::user_cache_dir(appname = "drugTargetInteractions")}) -
#' no Repurposing Hub data is bundled with or downloaded by the package
#' itself until this is called explicitly.
#'
#' If this fails with a certificate verification error, see the note in
#' this file's header comment (\code{broadRepurposingHubAccess.R}):
#' \code{repo-hub.broadinstitute.org} has occasionally served an
#' incomplete TLS certificate chain, a server-side issue unrelated to
#' this function.
#'
#' @param rerun logical(1); if \code{TRUE} (default), check for updates
#'   and (re)download as needed; if \code{FALSE}, use whatever is
#'   already in the local cache without checking upstream.
#' @param config list as returned by \code{genConfig()}; unused here but
#'   accepted for consistency with the rest of the package's API.
#' @return A named list of local file paths: \code{drug}, \code{sample}.
#' @examples
#' \donttest{
#'   paths <- downloadBroadRepurposingHub()
#'   paths
#' }
#' @seealso \code{\link{buildBroadRepurposingHubDb}}
#' @export
downloadBroadRepurposingHub <- function(rerun = TRUE, config = genConfig()) {
    ep    <- .brhEndpoints()
    paths <- vector("list", length(ep))
    names(paths) <- names(ep)
    for (nm in names(ep)) {
        fname <- basename(ep[[nm]])
        if (rerun) {
            paths[[nm]] <- .downloadFile(ep[[nm]], fname)
        } else {
            paths[[nm]] <- .getCacheFile(fname)
        }
    }
    paths
}


## ---------------------------------------------------------------------
## Explode multi-target drug rows into one row per (pert_iname, gene)
## ---------------------------------------------------------------------

#' Explode a drug table's " | "-separated \code{target} column into one
#' row per (\code{pert_iname}, gene) pair. Drugs with no listed target
#' (\code{NA}) keep exactly one row with \code{target_gene = NA}, so they
#' remain reachable via drug-name lookup even though they can never
#' surface via a gene-side query.
#' @keywords internal
.brhExplodeTargets <- function(drug) {
    genes <- strsplit(drug$target, " | ", fixed = TRUE)
    n <- lengths(genes)
    n[is.na(drug$target)] <- 1L
    genes[is.na(drug$target)] <- NA_character_
    data.frame(
        target_gene    = unlist(genes, use.names = FALSE),
        pert_iname     = rep(drug$pert_iname, n),
        clinical_phase = rep(drug$clinical_phase, n),
        moa            = rep(drug$moa, n),
        disease_area   = rep(drug$disease_area, n),
        indication     = rep(drug$indication, n),
        stringsAsFactors = FALSE
    )
}


## ---------------------------------------------------------------------
## Build (or reuse) a local SQLite, versioned from the raw files' own
## embedded "!File_date" - mirrors buildTtdDb()'s pattern of caching one
## ready-to-query SQLite file via BiocFileCache.
## ---------------------------------------------------------------------

#' Build (or fetch a cached) local SQLite database of Repurposing Hub
#' drug-target annotations
#'
#' Downloads the Repurposing Hub flat files (see
#' \code{\link{downloadBroadRepurposingHub}}), explodes the drug table's
#' multi-gene \code{target} column (one row per gene -
#' \code{\link{.brhExplodeTargets}}), left-joins in the sample-level
#' annotation (vendor/purity/structure), and writes a single
#' denormalized, indexed \code{broad_interactions} table to a local
#' SQLite file - the Repurposing Hub's target-drug relationship is as
#' flat as TTD's, so the same one-table shape applies (see
#' \code{\link{buildTtdDb}}). The file is cached via BiocFileCache under
#' a name that includes the source files' own \code{"!File_date"} (e.g.
#' \code{broad_repurposing_20250818.db}), so rebuilding is a no-op until
#' the Hub actually republishes. As with \code{\link{buildTtdDb}}, the
#' package itself never ships or redistributes this file - it is built
#' into the caller's own local cache the first time this is run.
#'
#' @param rerun logical(1); passed to
#'   \code{\link{downloadBroadRepurposingHub}}, and also controls whether
#'   an existing cached SQLite for the current version is reused
#'   (\code{FALSE}, the default here) or rebuilt (\code{TRUE}).
#'   Deliberately defaults to \code{FALSE} for the same reason
#'   \code{\link{buildTtdDb}} does: a bare \code{buildBroadRepurposingHubDb()}
#'   call is routine in examples/tests/vignette chunks, and the Hub
#'   updates rarely enough that defaulting to \code{TRUE} would just
#'   accumulate redundant cached copies with no benefit. Pass
#'   \code{rerun = TRUE} explicitly to force a fresh check for an update.
#' @param config list as returned by \code{genConfig()}.
#' @return character(1) local file path to the SQLite database.
#' @examples
#' \donttest{
#'   dbPath <- buildBroadRepurposingHubDb()
#'   dbPath
#' }
#' @seealso \code{\link{downloadBroadRepurposingHub}}, \code{\link{broadRepurposingHubAnnot}}
#' @export
buildBroadRepurposingHubDb <- function(rerun = FALSE, config = genConfig()) {
    paths <- downloadBroadRepurposingHub(rerun = rerun, config = config)
    drug  <- .brhReadTable(paths$drug)

    fileDate <- attr(drug, "fileDate")
    version <- if (!is.na(fileDate)) {
        d <- as.Date(fileDate, format = "%m/%d/%Y")
        if (is.na(d)) format(Sys.Date(), "%Y%m%d") else format(d, "%Y%m%d")
    } else {
        format(Sys.Date(), "%Y%m%d")
    }
    dbName <- paste0("broad_repurposing_", version, ".db")

    if (!rerun) {
        existing <- tryCatch(.getCacheFile(dbName), error = function(e) NA_character_)
        if (length(existing) > 0 && !is.na(existing)) return(existing)
    }

    sample <- .brhReadTable(paths$sample)
    sample <- sample[!duplicated(sample), ]

    exploded <- .brhExplodeTargets(drug)
    interactions <- merge(exploded, sample, by = "pert_iname", all.x = TRUE)
    interactions <- interactions[, c("target_gene", "pert_iname", "clinical_phase",
                                     "moa", "disease_area", "indication", "broad_id",
                                     "qc_incompatible", "purity", "vendor", "catalog_no",
                                     "vendor_name", "expected_mass", "smiles", "InChIKey",
                                     "pubchem_cid", "deprecated_broad_id")]

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "broad_interactions", interactions, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_brh_target   ON broad_interactions (target_gene)")
    dbExecute(con, "CREATE INDEX idx_brh_pert     ON broad_interactions (pert_iname)")
    dbExecute(con, "CREATE INDEX idx_brh_broad_id ON broad_interactions (broad_id)")
    dbDisconnect(con)

    bfc <- .getCache()
    rid <- names(bfcadd(bfc, dbName, tmpDb, action = "copy"))
    file.remove(tmpDb)
    bfcrpath(bfc, rids = rid)
}


## ---------------------------------------------------------------------
## Bidirectional query, mirroring ttdTargetAnnot()'s queryBy convention
## ---------------------------------------------------------------------

#' Query Broad Repurposing Hub drug-target annotations bidirectionally
#'
#' Queries the local Repurposing Hub SQLite built by
#' \code{\link{buildBroadRepurposingHubDb}}, using the same
#' \code{queryBy = list(molType, idType, ids)} convention as
#' \code{\link{ttdTargetAnnot}} - including its \code{QueryIDs} column:
#' every row of the result is tagged with the original query token it
#' matched, and query IDs that returned no rows still appear as a single
#' row with all other fields \code{NA}.
#' \itemize{
#'   \item \code{molType = "protein"} (or \code{"gene"}), \code{idType =
#'     "symbol"} -> target -> drug. The Repurposing Hub only exposes
#'     gene symbols as targets (no accession/ID system of its own).
#'   \item \code{molType = "cmp"}, \code{idType} one of \code{"name"}
#'     (\code{pert_iname}) or \code{"broad_id"} (a specific physical
#'     sample/batch ID) -> drug -> target.
#' }
#' \code{"symbol"}/\code{"name"} lookups are case-insensitive (the
#' Repurposing Hub's own \code{pert_iname} values are lower-case, unlike
#' most other sources in this package); \code{"broad_id"} is exact-match.
#' A drug with no listed target still has one row with
#' \code{target_gene = NA} (see \code{\link{.brhExplodeTargets}}), so it
#' remains reachable by name/broad_id even though it can never surface
#' via a gene-side query.
#'
#' @param queryBy list with components \code{molType}, \code{idType},
#'   \code{ids} (character vector).
#' @param brhDbPath character(1) path to the Repurposing Hub SQLite, e.g.
#'   from \code{\link{buildBroadRepurposingHubDb}}.
#' @param fields \code{"core"} (default) or \code{"all"} - both return
#'   every \code{broad_interactions} column, since (unlike the REST-backed
#'   sources) there is no larger raw payload to opt into - or a character
#'   vector of column names to keep (\code{QueryIDs} is always retained).
#'   See \code{\link{listDrugTargetFields}}.
#' @return A \code{data.frame} with columns \code{QueryIDs},
#'   \code{target_gene}, \code{pert_iname}, \code{clinical_phase},
#'   \code{moa}, \code{disease_area}, \code{indication}, \code{broad_id},
#'   \code{qc_incompatible}, \code{purity}, \code{vendor},
#'   \code{catalog_no}, \code{vendor_name}, \code{expected_mass},
#'   \code{smiles}, \code{InChIKey}, \code{pubchem_cid},
#'   \code{deprecated_broad_id} (with the default \code{fields = "core"}),
#'   or a subset when \code{fields} requests specific columns.
#' @examples
#' \donttest{
#'   dbPath <- buildBroadRepurposingHubDb()
#'   broadRepurposingHubAnnot(list(molType = "protein", idType = "symbol",
#'                                 ids = c("FGFR1", "IL1B")), dbPath)
#'   broadRepurposingHubAnnot(list(molType = "cmp", idType = "name",
#'                                 ids = "pemigatinib"), dbPath)
#' }
#' @seealso \code{\link{buildBroadRepurposingHubDb}}, \code{\link{ttdTargetAnnot}},
#'   \code{\link{listDrugTargetFields}}
#' @export
broadRepurposingHubAnnot <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                                     brhDbPath, fields = "core") {
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
               symbol = "target_gene",
               stop("idType for molType='protein'/'gene' must be 'symbol'"))
    } else if (queryBy$molType == "cmp") {
        switch(queryBy$idType,
               name     = "pert_iname",
               broad_id = "broad_id",
               stop("idType for molType='cmp' must be one of: ",
                    "'name', 'broad_id'"))
    } else {
        stop("molType must be 'protein'/'gene' or 'cmp'")
    }

    caseInsensitive <- queryBy$idType %in% c("symbol", "name")
    ids <- if (caseInsensitive) toupper(queryBy$ids) else queryBy$ids
    idvec <- paste0("('", paste(gsub("'", "''", ids, fixed = TRUE), collapse = "', '"), "')")
    colExpr <- if (caseInsensitive) paste0("UPPER(", col, ")") else col

    con <- dbConnect(SQLite(), brhDbPath)
    on.exit(dbDisconnect(con))
    query <- paste0("SELECT * FROM broad_interactions WHERE ", colExpr, " IN ", idvec)
    resultDF <- dbGetQuery(con, query)

    ## Tag every row with the original query token it matched, and NA-pad
    ## any query ID that matched nothing - mirrors ttdTargetAnnot()'s
    ## QueryIDs convention exactly (same Inf-index trick: an unmatched ID
    ## gets rowid Inf, and indexing a data.frame with Inf yields a row of
    ## NAs).
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
    .dtiSelectFields(out, fields, .dtiBroadAllCols)
}

#' Documented column list for \code{listDrugTargetFields("broad")}
#'
#' Unlike the REST-backed sources' \code{fields = "all"} (ChEMBL,
#' PubChem, DGIdb, Open Targets), \code{broad_interactions} is a single
#' flat local SQLite table built entirely by
#' \code{\link{buildBroadRepurposingHubDb}} (see there), so this is an
#' exact list, not a best-effort one, and \code{fields = "core"} and
#' \code{fields = "all"} are equivalent for
#' \code{\link{broadRepurposingHubAnnot}}.
#' @keywords internal
.dtiBroadAllCols <- c("QueryIDs", "target_gene", "pert_iname", "clinical_phase",
                      "moa", "disease_area", "indication", "broad_id",
                      "qc_incompatible", "purity", "vendor", "catalog_no",
                      "vendor_name", "expected_mass", "smiles", "InChIKey",
                      "pubchem_cid", "deprecated_broad_id")
