## =====================================================================
##  unichemAccess.R
##  Compound-ID cross-referencing for the drugTargetInteractions
##  Bioconductor package, via a local SQLite built from UniChem's own
##  bulk table dumps. Mirrors downloadChemblDb()/buildTtdDb()'s
##  local-SQL pattern (download once, build a local database, query it
##  repeatedly) rather than apiAccess.R's/idTranslation.R's live-API
##  pattern, for the same reason TTD does: UniChem's live per-compound
##  REST API dropped batch-query support in its 2.0 relaunch and (live-
##  tested 2026-07-17) is intermittently unreliable even for single
##  lookups spaced several seconds apart - fine for the one-off,
##  low-volume lookups already covered by other functions, not for
##  building a comprehensive cross-reference table.
##
##  Design goal: a single, source-agnostic cross-reference table, not a
##  hand-coded translation per database pair. UniChem publishes two
##  different bulk-download mechanisms for different jobs (confirmed
##  2026-07-17, both regenerated together from the same underlying
##  database, neither superseding the other):
##   - `wholeSourceMapping/` - per-source-PAIR files (e.g. src1src22.txt.gz
##     = ChEMBL<->PubChem), meant for one-off bulk conversion between two
##     *specific* named databases. Covering every source generically this
##     way would need up to C(24,2) = 253 separate downloads.
##   - `table_dumps/reference.tsv.gz` - UniChem's own internal reference
##     table, UCI-centric: one row per (compound, source), sharing a
##     common UCI (UniChem Identifier) across every source a given
##     structure appears in. This is the "build your own local index"
##     format, and is what this file uses - a single self-join on UCI
##     gives full any-source-to-any-source translation from one file.
##  That is what makes getUnichemMapping()'s `from`/`to` arguments accept
##  *any* source name UniChem carries - including ones this package has
##  no dedicated wrapper for (e.g. CCDC, live-confirmed working via this
##  exact mechanism on 2026-07-17) - entirely generically, with zero
##  source-specific code. Re-running buildUnichemDb() against a future
##  UniChem refresh picks up newly added sources (or a previously
##  dropped one, e.g. LINCS, was a UniChem source in the past but is not
##  among the current 24 as of 2026-07-17) automatically.
##
##  Only reference.tsv.gz (~1.5GB compressed, ~166M rows across 24
##  sources) and source.tsv.gz (~3KB, source-ID/name lookup) are used;
##  UniChem's third table dump, structure.tsv.gz (~13.4GB, actual
##  chemical structures/InChI), is deliberately NOT downloaded - it is
##  not needed for ID-to-ID translation.
##
##  Scope: rows are kept only for UCIs that (a) appear under >=
##  `minSources` distinct sources (default 2 - a UCI present in exactly
##  one source has nothing to translate to or from), AND (b) include at
##  least one source from `anchorSources` (default: ChEMBL, DrugBank,
##  ChEBI, Guide to Pharmacology, DrugCentral, BindingDB, ClinicalTrials
##  - see .unichemAnchorSources). Condition (a) alone is not enough:
##  live-tested 2026-07-17, building the unfiltered (a)-only version
##  first kept 79.5M of 172M rows, but 87% of those came from
##  pubchem/surechembl/molport clusters with zero connection to anything
##  drug-relevant (PubChem ingests both as substance providers, so they
##  cross-reference each other heavily on their own). Condition (b)
##  brings that down to ~13.1M rows / ~730MB while still keeping every
##  PubChem (or any other source's) entry that *does* cross-reference
##  something drug-relevant - see .unichemFilterAndIndex() for the full
##  reasoning and .unichemAnchorSources for the exact list. A PubChem-only
##  compound with real target/bioactivity data (no drug-relevant
##  cross-reference at all) is still fully queryable via its native
##  PubChem CID through getPubchemDrugs()/getPubchemTargets() in
##  apiAccess.R - that is compound-target *discovery*, a different job
##  from this file's *translation* role, and does not depend on this
##  table at all.
##
##  Distribution posture (interim, 2026-07-17): the built SQLite is
##  cached locally via the package's existing BiocFileCache
##  infrastructure, same as downloadChemblDb()/buildTtdDb(). Until this
##  is registered as a proper Bioconductor ExperimentHub resource
##  (targeted before the October Bioc release), the pre-built file is
##  distributed from the same S3 bucket used for the legacy 3-file
##  UniChem mirror (see downloadUniChem() in
##  drugTargetAnnotations_Fct.R) - manual upload step, not automated by
##  this file.
##
##  Suggested DESCRIPTION addition: none (uses only RSQLite/BiocFileCache,
##  already Imports).
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for UniChem's bulk table dumps
#' @keywords internal
.unichemEndpoints <- function() {
    list(
        base      = "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/table_dumps",
        reference = "reference.tsv.gz",
        source    = "source.tsv.gz"
    )
}


## ---------------------------------------------------------------------
## Download raw UniChem table dumps (cached via the package's existing
## BiocFileCache helpers)
## ---------------------------------------------------------------------

#' Download the UniChem table dumps needed for compound cross-referencing
#'
#' Downloads (or reuses previously cached copies of) UniChem's
#' \code{reference} (UCI/source/compound-ID triples, ~1.5GB) and
#' \code{source} (source-ID/name lookup, ~3KB) table dumps via the
#' package's existing \code{.downloadFile()}/BiocFileCache
#' infrastructure - the same mechanism \code{downloadChemblDb()} and
#' \code{downloadTTD()} use. \code{structure.tsv.gz} (actual chemical
#' structures, ~13.4GB) is deliberately not downloaded - see the file
#' header.
#'
#' @param rerun logical(1); if \code{TRUE} (default), check for updates
#'   and (re)download as needed; if \code{FALSE}, use whatever is
#'   already in the local cache without checking upstream.
#' @param config list as returned by \code{genConfig()}; unused here but
#'   accepted for consistency with the rest of the package's API.
#' @return A named list of local file paths: \code{reference},
#'   \code{source}.
#' @examples
#' \donttest{
#'   paths <- downloadUnichemTables()
#'   paths
#' }
#' @seealso \code{\link{buildUnichemDb}}
#' @export
downloadUnichemTables <- function(rerun = TRUE, config = genConfig()) {
    ep    <- .unichemEndpoints()
    files <- ep[c("reference", "source")]
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
## Build (or reuse) a local SQLite cross-reference table
## ---------------------------------------------------------------------

#' Chunked-load UniChem's reference table dump into a staging SQLite table
#'
#' Streams \code{reference.tsv.gz} through a gzip connection in batches
#' of \code{chunkSize} rows (never materialising the full ~166M-row file
#' as an in-memory data.frame) and appends each batch to a
#' \code{ref_raw} table. \code{SRC_COMPOUND_ID} is read as character
#' throughout - several sources (e.g. SwissLipids' \code{"SLM:..."} IDs)
#' are not purely numeric.
#'
#' @param refPath character(1) path to the downloaded \code{reference.tsv.gz}.
#' @param con an open \code{RSQLite} connection.
#' @param chunkSize integer(1) rows read per batch (default 500,000;
#'   live-tested 2026-07-17 at ~4s/2.3M-row-slice, i.e. ~5 minutes
#'   extrapolated for the full file).
#' @param verbose logical(1); if TRUE, message progress per chunk.
#' @return integer(1) total rows loaded (invisibly).
#' @keywords internal
.unichemLoadReference <- function(refPath, con, chunkSize = 500000L, verbose = FALSE) {
    f <- gzfile(refPath, "rt")
    on.exit(close(f))
    invisible(readLines(f, n = 1L))  ## header
    totalRows <- 0L
    repeat {
        chunk <- tryCatch(
            utils::read.delim(f, header = FALSE, nrows = chunkSize, sep = "\t",
                              col.names = c("uci", "src_id", "compound_id", "assignment"),
                              colClasses = c("integer", "integer", "character", "integer"),
                              quote = "", comment.char = ""),
            error = function(e) NULL)
        if (is.null(chunk) || nrow(chunk) == 0L) break
        dbWriteTable(con, "ref_raw", chunk[, c("uci", "src_id", "compound_id")],
                          append = TRUE)
        totalRows <- totalRows + nrow(chunk)
        if (verbose) message("UniChem reference: loaded ", totalRows, " rows so far")
        if (nrow(chunk) < chunkSize) break
    }
    invisible(totalRows)
}

#' Default anchor sources for \code{\link{buildUnichemDb}}
#'
#' Drug/bioactivity-relevant UniChem sources. See
#' \code{\link{.unichemFilterAndIndex}} for why plain "\code{>=2
#' sources}" is not a strict enough filter on its own - live-tested
#' 2026-07-17: it retained 79.5M of 172M rows, 87% of which came from
#' just pubchem/surechembl/molport clusters with no connection to any
#' of these anchor sources at all.
#' @keywords internal
.unichemAnchorSources <- c("chembl", "drugbank", "chebi", "gtopdb",
                           "drugcentral", "bindingdb", "clinicaltrials")

#' Filter the staging table to multi-source, anchor-relevant UCIs and
#' build indexes
#'
#' Keeps only rows belonging to a UCI (UniChem structure cluster) that
#' (a) appears under at least \code{minSources} distinct sources, and
#' (b) includes at least one row from \code{anchorSources} (resolved via
#' the already-loaded \code{unichem_sources} table). Condition (a) alone
#' is not enough to keep this table drug-target-relevant: PubChem cross-
#' references SureChEMBL (patent-mined structures) and MolPort (a
#' chemical vendor's catalog) heavily, so ">=2 sources" alone let through
#' huge PubChem<->SureChEMBL<->MolPort clusters with zero connection to
#' ChEMBL, DrugBank, or any other source this package's drug-target work
#' actually cares about (confirmed live 2026-07-17 by building the
#' unfiltered version first: 79.5M of 172M rows survived ">=2 sources",
#' but 87% of those rows came from pubchem/surechembl/molport alone).
#' Condition (b) fixes this without discarding PubChem entries wholesale
#' - a PubChem compound is kept exactly when it also cross-references
#' something drug-relevant, which is the actually-useful case for a
#' *translation* table (see the file header for why PubChem-only compound-
#' target data doesn't need this table at all).
#'
#' Applied as a join against a temporary \code{anchor_ucis} table rather
#' than a large SQL \code{IN (...)} subquery - live-tested 2026-07-17:
#' the subquery form did not complete in 2 minutes against 79.5M rows;
#' the join form did.
#'
#' @param con an open \code{RSQLite} connection with a populated
#'   \code{ref_raw} table and (if \code{anchorSources} is not
#'   \code{NULL}) an already-loaded \code{unichem_sources} table.
#' @param minSources integer(1) minimum distinct sources per UCI to keep
#'   (default 2).
#' @param anchorSources character vector of UniChem source names (see
#'   \code{\link{.unichemAnchorSources}} for the default), or \code{NULL}
#'   to disable the anchor requirement and keep condition (a) only.
#' @param verbose logical(1); if TRUE, message progress.
#' @return invisible(NULL).
#' @keywords internal
.unichemFilterAndIndex <- function(con, minSources = 2L,
                                   anchorSources = .unichemAnchorSources,
                                   verbose = FALSE) {
    if (verbose) message("UniChem: indexing staging table")
    dbExecute(con, "CREATE INDEX idx_raw_uci ON ref_raw (uci)")

    if (verbose) message("UniChem: filtering to UCIs with >= ", minSources, " sources")
    dbExecute(con, sprintf(
        "CREATE TABLE unichem_xref AS
         SELECT uci, src_id, compound_id FROM ref_raw
         WHERE uci IN (SELECT uci FROM ref_raw GROUP BY uci HAVING COUNT(DISTINCT src_id) >= %d)",
        minSources))
    dbExecute(con, "DROP TABLE ref_raw")
    dbExecute(con, "CREATE INDEX idx_xref_uci ON unichem_xref (uci)")

    if (!is.null(anchorSources)) {
        if (verbose) message("UniChem: restricting to UCIs referencing >= 1 of: ",
                             paste(anchorSources, collapse = ", "))
        anchorIds <- vapply(anchorSources, .unichemResolveSourceId, integer(1), con = con)
        dbExecute(con, sprintf(
            "CREATE TABLE anchor_ucis AS
             SELECT DISTINCT uci FROM unichem_xref WHERE src_id IN (%s)",
            paste(anchorIds, collapse = ",")))
        dbExecute(con, "CREATE INDEX idx_anchor_uci ON anchor_ucis (uci)")
        dbExecute(con,
            "CREATE TABLE unichem_xref2 AS
             SELECT x.uci, x.src_id, x.compound_id FROM unichem_xref x
             JOIN anchor_ucis a ON x.uci = a.uci")
        dbExecute(con, "DROP TABLE unichem_xref")
        dbExecute(con, "DROP TABLE anchor_ucis")
        dbExecute(con, "ALTER TABLE unichem_xref2 RENAME TO unichem_xref")
        dbExecute(con, "CREATE INDEX idx_xref_uci ON unichem_xref (uci)")
    }

    if (verbose) message("UniChem: building final indexes")
    dbExecute(con, "CREATE INDEX idx_xref_src_cmp ON unichem_xref (src_id, compound_id)")
    dbExecute(con, "VACUUM")
    invisible(NULL)
}

#' Build (or fetch a cached) local SQLite of UniChem compound cross-references
#'
#' Downloads UniChem's bulk table dumps (see
#' \code{\link{downloadUnichemTables}}), loads \code{reference.tsv.gz} in
#' chunks (never all in memory at once), filters to compounds with a
#' cross-reference to at least \code{minSources} distinct databases, and
#' writes two tables to a local SQLite file: \code{unichem_xref(uci,
#' src_id, compound_id)} and \code{unichem_sources(src_id, name,
#' name_long, description)} (a direct copy of \code{source.tsv.gz}, used
#' by \code{\link{getUnichemMapping}} to resolve source-name arguments).
#' The file is cached via BiocFileCache under a name that includes the
#' build date (UniChem's bulk dumps carry no release-version string the
#' way TTD's flat files do), so rebuilding on the same day is a no-op
#' unless \code{rerun = TRUE}.
#'
#' This is a substantially heavier operation than the package's other
#' \code{build*Db()}/\code{download*Db()} functions - live-tested
#' 2026-07-17: ~1.5GB download, ~166M source rows to load/filter/index.
#' Expect this to take real time (minutes, not seconds) and several GB
#' of temporary disk space for the unfiltered staging table before
#' filtering shrinks it down.
#'
#' @param rerun logical(1); passed to \code{\link{downloadUnichemTables}},
#'   and also controls whether an existing cached SQLite for today's
#'   build-date is reused (\code{FALSE}) or rebuilt (\code{TRUE}).
#' @param config list as returned by \code{genConfig()}.
#' @param minSources integer(1) minimum distinct sources per UCI to keep
#'   (default 2) - see \code{\link{.unichemFilterAndIndex}}.
#' @param anchorSources character vector of UniChem source names a kept
#'   UCI must reference at least one of, or \code{NULL} to disable this
#'   requirement (default: \code{\link{.unichemAnchorSources}} - ChEMBL,
#'   DrugBank, ChEBI, Guide to Pharmacology, DrugCentral, BindingDB,
#'   ClinicalTrials). See \code{\link{.unichemFilterAndIndex}} for why
#'   this matters: \code{minSources} alone is not enough to keep the
#'   table drug-target-relevant.
#' @param chunkSize integer(1) rows read per batch while loading
#'   \code{reference.tsv.gz} (default 500,000).
#' @param verbose logical(1); if TRUE, message progress throughout.
#' @return character(1) local file path to the SQLite database.
#' @examples
#' \donttest{
#'   dbPath <- buildUnichemDb(verbose = TRUE)
#'   dbPath
#' }
#' @seealso \code{\link{downloadUnichemTables}}, \code{\link{getUnichemMapping}}
#' @export
buildUnichemDb <- function(rerun = TRUE, config = genConfig(), minSources = 2L,
                           anchorSources = .unichemAnchorSources,
                           chunkSize = 500000L, verbose = FALSE) {
    dbName <- paste0("unichem_", format(Sys.Date(), "%Y%m%d"), ".db")

    if (!rerun) {
        existing <- tryCatch(.getCacheFile(dbName), error = function(e) NA_character_)
        if (length(existing) > 0 && !is.na(existing)) return(existing)
    }

    paths <- downloadUnichemTables(rerun = rerun, config = config)

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    on.exit(dbDisconnect(con), add = TRUE)
    dbExecute(con, "PRAGMA journal_mode=OFF")
    dbExecute(con, "PRAGMA synchronous=OFF")

    ## Sources are loaded first: anchor-source resolution in
    ## .unichemFilterAndIndex() needs unichem_sources to already exist.
    sources <- utils::read.delim(gzfile(paths$source), header = TRUE, sep = "\t",
                                 col.names = c("src_id", "name", "name_long", "description"),
                                 colClasses = "character", quote = "", comment.char = "")
    sources$src_id <- as.integer(sources$src_id)
    dbWriteTable(con, "unichem_sources", sources, overwrite = TRUE)

    if (verbose) message("UniChem: loading reference table (this takes a while)")
    .unichemLoadReference(paths$reference, con, chunkSize = chunkSize, verbose = verbose)
    .unichemFilterAndIndex(con, minSources = minSources, anchorSources = anchorSources,
                           verbose = verbose)

    dbDisconnect(con)
    on.exit(NULL)  ## already disconnected above

    bfc <- .getCache()
    rid <- names(bfcadd(bfc, dbName, tmpDb, action = "copy"))
    file.remove(tmpDb)
    bfcrpath(bfc, rids = rid)
}


## ---------------------------------------------------------------------
## Generic query, source-agnostic by construction
## ---------------------------------------------------------------------

#' Resolve a UniChem \code{from}/\code{to} argument to a numeric source ID
#'
#' Case-insensitive on the source \code{name}: UniChem's own
#' \code{source.tsv.gz} is not consistently cased (most names are
#' lowercase, e.g. \code{"chembl"}, but at least one, \code{"CCDC"}, is
#' not - live-confirmed 2026-07-17), so matching is done via
#' \code{LOWER(name)} rather than assuming any particular casing.
#' @keywords internal
.unichemResolveSourceId <- function(x, con) {
    if (grepl("^[0-9]+$", x)) return(as.integer(x))
    hit <- dbGetQuery(con, "SELECT src_id FROM unichem_sources WHERE LOWER(name) = ?",
                           params = list(tolower(x)))
    if (nrow(hit) == 0L)
        stop("Unrecognised UniChem source '", x, "'. See unichem_sources in the ",
             "database (or https://www.ebi.ac.uk/unichem/api/v1/sources/) for the ",
             "current list of valid names.")
    hit$src_id[1]
}

#' Translate compound identifiers between databases via a local UniChem SQLite
#'
#' Genuinely source-agnostic: \code{from}/\code{to} accept any source
#' name (or raw numeric source ID) present in the \code{unichem_sources}
#' table of the database built by \code{\link{buildUnichemDb}} - not a
#' fixed enum hand-coded per database. A source this package has no
#' dedicated accessor for (e.g. \code{"gtopdb"}, \code{"drugcentral"},
#' or any source added to a future UniChem refresh) works identically to
#' \code{"chembl"} or \code{"pubchem"} the moment it is present in the
#' underlying data - no code changes needed. Implemented as a single SQL
#' self-join on the shared UniChem structure cluster ID (UCI), which is
#' what makes this generic: the query never references a specific source
#' by name in its logic, only in the two parameter values.
#'
#' @param ids character vector of source-database identifiers.
#' @param from character(1) source database name (e.g. \code{"chembl"},
#'   \code{"pubchem"}, \code{"drugbank"}, \code{"chebi"}) or numeric
#'   source ID, matching \code{unichem_sources$name}/\code{src_id}.
#' @param to character(1) target database name or numeric source ID.
#' @param dbPath character(1) path to the UniChem SQLite, e.g. from
#'   \code{\link{buildUnichemDb}}.
#' @return A \code{data.frame} with columns \code{From} and \code{To}.
#'   IDs with no cross-reference to \code{to} contribute no row (the
#'   underlying table only carries compounds with >= 2 source
#'   cross-references in the first place - see
#'   \code{\link{buildUnichemDb}} - so an unmapped ID may mean either no
#'   \code{to}-side entry exists, or the compound had no cross-reference
#'   to *any* other source and was excluded from the table entirely).
#' @examples
#' \donttest{
#'   dbPath <- buildUnichemDb()
#'   getUnichemMapping("CHEMBL25", from = "chembl", to = "pubchem", dbPath)
#'   getUnichemMapping("CHEMBL25", from = "chembl", to = "drugbank", dbPath)
#' }
#' @seealso \code{\link{buildUnichemDb}}, \code{\link{getUniprotMapping}}
#' @export
getUnichemMapping <- function(ids, from, to, dbPath) {
    stopifnot(is.character(ids), length(ids) >= 1L,
             is.character(from) || is.numeric(from), length(from) == 1L,
             is.character(to) || is.numeric(to), length(to) == 1L)
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))

    fromId <- .unichemResolveSourceId(as.character(from), con)
    toId   <- .unichemResolveSourceId(as.character(to), con)

    idvec <- paste0("('", paste(gsub("'", "''", unique(ids), fixed = TRUE),
                                collapse = "', '"), "')")
    query <- sprintf(
        "SELECT a.compound_id AS 'From', b.compound_id AS 'To'
         FROM unichem_xref a JOIN unichem_xref b ON a.uci = b.uci
         WHERE a.src_id = %d AND a.compound_id IN %s AND b.src_id = %d",
        fromId, idvec, toId)
    dbGetQuery(con, query)
}
