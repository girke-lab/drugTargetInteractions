## =====================================================================
##  genomeWideAnnot.R
##  Genome-wide, HGNC-anchored master drug-target annotation table for
##  the drugTargetInteractions Bioconductor package.
##
##  Motivation: queryDrugTargets()/combineDrugTargets() (drugTargetMeta.R)
##  are designed for on-demand queries of a handful to a few thousand
##  IDs. Answering genome-wide questions (e.g. "do disease-associated
##  gene variants have known drugs annotated?", "do LINCS perturbation
##  signature hits have known targets/MOAs?") needs a periodically
##  rebuilt, shareable master table covering all ~19,200 human
##  protein-coding genes instead - a fundamentally different scale, with
##  its own concerns (checkpointing/resumability, a canonical gene list,
##  symbol-drift reconciliation) that this file adds as first-class
##  pieces rather than stretching the on-demand dispatcher to cover both.
##
##  A single target-centric run across the 4 curated annotation sources
##  (ChEMBL, DGIdb, Open Targets, TTD) already answers the reverse
##  drug -> target question too: every row already carries whichever
##  drug matched that gene, so filtering the assembled table by drug
##  identifier gives that direction for free. ChEMBL's bioassay function
##  and PubChem are deliberately NOT part of this build - see
##  getChemblBioassay()/getPubchemDrugTarget() (apiAccess.R) for the
##  separate raw-measurement track, which is a different kind of job at
##  genome scale (bulk downloads, not this live-API loop).
## =====================================================================


## ---------------------------------------------------------------------
## HGNC gene table (the gene-centric anchor)
## ---------------------------------------------------------------------
## HGNC (https://www.genenames.org/) publishes one row per approved
## human gene with the current official symbol, all previous/alias
## symbols, and cross-references to Ensembl/UniProt/Entrez in a single
## file - chosen over EnsDb.Hsapiens.v86 or a live UniProt query
## specifically because it is the single authoritative source for
## symbol *and* prev_symbol/alias_symbol together, which the DGIdb/TTD
## symbol-drift reconciliation below depends on. Files are hosted on a
## Google Cloud Storage bucket (the older EBI FTP path is defunct).

#' Pinned HGNC quarterly archive filename
#'
#' HGNC archives quarterly snapshots at irregular dates (not simply
#' YYYY-MM-01) - confirmed live 2026-07-18 via the bucket's listing API:
#' \code{hgnc_complete_set_2026-01-06.txt}, \code{_2026-04-01.txt},
#' \code{_2026-04-07.txt}, \code{_2026-07-03.txt}, \code{_2026-07-07.txt}
#' all exist; \code{_2026-07-01.txt} does not. A fixed, confirmed-working
#' filename is pinned here as the package default so vignette/test
#' output and downstream analyses stay reproducible across time (see
#' \code{current = TRUE} in \code{\link{downloadHgncTable}} for the
#' always-latest rolling file instead). Update this string periodically
#' as a maintenance task - check availability first via
#' \url{https://storage.googleapis.com/public-download-files?prefix=hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_}.
#' @keywords internal
#' @noRd
.hgncPinnedQuarterlyFile <- "hgnc_complete_set_2026-04-01.txt"

#' Download (and cache) an HGNC complete-gene-set TSV snapshot
#'
#' @param archiveFile character(1) or \code{NULL}; a specific HGNC
#'   archive filename (quarterly or monthly, see
#'   \url{https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/})
#'   to reproduce a particular snapshot. \code{NULL} (default) uses the
#'   package's pinned quarterly snapshot (see
#'   \code{.hgncPinnedQuarterlyFile}) - fixed across package
#'   versions so results stay reproducible; the pin is updated as a
#'   package maintenance task, not a user setting.
#' @param current logical(1); if \code{TRUE}, ignore \code{archiveFile}
#'   and download HGNC's always-latest rolling file instead (updated
#'   Tuesdays/Fridays) - convenient, but re-running the same code later
#'   can yield different results.
#' @param rerun logical(1); if \code{TRUE} (default), check for a newer
#'   remote version even if already cached (only matters for
#'   \code{current = TRUE}; pinned/archive snapshots never change once
#'   published). If \code{FALSE}, reuse whatever is cached without
#'   checking.
#' @param config list from \code{\link{genConfig}}.
#' @param verbose logical(1).
#' @return character(1) local file path to the cached TSV.
#' @examples
#' \donttest{
#'   hgncPath <- downloadHgncTable()
#' }
#' @seealso \code{\link{getHgncGeneTable}}
#' @export
downloadHgncTable <- function(archiveFile = NULL, current = FALSE,
                              rerun = TRUE, config = genConfig(),
                              verbose = FALSE) {
    base <- "https://storage.googleapis.com/public-download-files/hgnc"
    if (isTRUE(current)) {
        name <- "hgnc_complete_set.txt"
        url <- paste0(base, "/tsv/tsv/", name)
    } else {
        name <- archiveFile %||% .hgncPinnedQuarterlyFile
        url <- paste0(base, "/archive/archive/quarterly/tsv/", name)
    }
    if (!rerun) {
        existing <- tryCatch(.getCacheFile(name), error = function(e) NA_character_)
        if (length(existing) > 0 && !is.na(existing)) return(existing)
    }
    .downloadFile(url, name, verbose = verbose)
}

#' Load the HGNC complete gene set as a gene-centric data.frame
#'
#' Reads the cached HGNC TSV (see \code{\link{downloadHgncTable}}),
#' optionally restricts to protein-coding genes (~19,200 of them), and
#' reshapes the pipe-separated multi-value columns
#' (\code{prev_symbol}/\code{alias_symbol}/\code{uniprot_ids}) into
#' list-columns - one row per approved HGNC gene report. Isoform-suffixed
#' UniProt accessions (e.g. \code{"P12345-2"}) are stripped defensively,
#' though HGNC's own cross-references are already canonical-only.
#'
#' @param proteinCodingOnly logical(1); if \code{TRUE} (default), keep
#'   only rows where \code{locus_group == "protein-coding gene"}.
#' @param archiveFile,current,rerun,config,verbose passed to
#'   \code{\link{downloadHgncTable}}; \code{rerun} defaults to
#'   \code{FALSE} here (unlike \code{downloadHgncTable()}'s own default)
#'   since this higher-level function is the one called routinely.
#' @return A \code{data.frame} with columns \code{hgnc_id}, \code{symbol},
#'   \code{prev_symbol} (list-column), \code{alias_symbol} (list-column),
#'   \code{entrez_id}, \code{ensembl_gene_id}, \code{uniprot_ids}
#'   (list-column), \code{locus_group}. \code{attr(., "hgncSource")}
#'   records the source filename for provenance.
#' @examples
#' \donttest{
#'   hgncTable <- getHgncGeneTable()
#'   nrow(hgncTable)
#'   hgncTable[hgncTable$symbol == "FGFR1", c("symbol", "uniprot_ids")]
#' }
#' @seealso \code{\link{downloadHgncTable}}, \code{\link{buildHgncSymbolMap}},
#'   \code{\link{buildGenomeWideDrugTargetTable}}
#' @export
getHgncGeneTable <- function(proteinCodingOnly = TRUE, archiveFile = NULL,
                             current = FALSE, rerun = FALSE,
                             config = genConfig(), verbose = FALSE) {
    path <- downloadHgncTable(archiveFile = archiveFile, current = current,
                              rerun = rerun, config = config, verbose = verbose)
    ## Default quote handling (not quote = "") - HGNC's TSV relies on
    ## standard CSV-style quoting to wrap fields whose own content needs
    ## it (confirmed live: some alias_symbol values are wrapped this way);
    ## disabling quote interpretation leaves literal quote characters in
    ## the data instead of stripping them.
    df <- read.delim(path, stringsAsFactors = FALSE,
                     na.strings = "", check.names = FALSE)
    if (isTRUE(proteinCodingOnly))
        df <- df[!is.na(df$locus_group) & df$locus_group == "protein-coding gene", ]

    keep <- c("hgnc_id", "symbol", "prev_symbol", "alias_symbol",
             "entrez_id", "ensembl_gene_id", "uniprot_ids", "locus_group")
    df <- df[, intersect(keep, names(df)), drop = FALSE]

    splitCol <- function(x) strsplit(ifelse(is.na(x), "", x), "\\|")
    if ("prev_symbol" %in% names(df)) df$prev_symbol <- splitCol(df$prev_symbol)
    if ("alias_symbol" %in% names(df)) df$alias_symbol <- splitCol(df$alias_symbol)
    if ("uniprot_ids" %in% names(df)) {
        df$uniprot_ids <- lapply(splitCol(df$uniprot_ids), function(acc) {
            acc <- acc[nzchar(acc)]
            unique(sub("-[0-9]+$", "", acc))
        })
    }
    rownames(df) <- NULL
    attr(df, "hgncSource") <- basename(path)
    df
}

#' Build a prev/alias-symbol -> current-symbol normalization map
#'
#' DGIdb and TTD (and potentially any other external gene list) echo
#' back whatever symbol was current when \emph{their} records were
#' curated; a nontrivial fraction are now outdated relative to HGNC.
#' This builds a lookup from every historical \code{prev_symbol}/
#' \code{alias_symbol} value to its gene's current approved
#' \code{symbol}, so incoming symbols from those sources can be routed
#' through it before joining back to an HGNC-anchored table.
#'
#' @param hgncTable data.frame as returned by \code{\link{getHgncGeneTable}}
#'   (must have \code{symbol}, \code{prev_symbol}, \code{alias_symbol}
#'   list-columns).
#' @param warn logical(1); if \code{TRUE} (default), raise a
#'   \code{warning()} reporting how many of the map's entries are
#'   ambiguous. That count describes the HGNC snapshot as a whole, not
#'   any particular set of symbols, so functions that build a map only to
#'   translate a handful of symbols pass \code{FALSE} and let
#'   \code{\link{normalizeGeneSymbols}} report the ambiguities those
#'   symbols actually run into.
#' @return A named character vector (old symbol -> current symbol).
#'   Old symbols that map to more than one current symbol (real cases
#'   exist, e.g. shared paralog aliases) keep only the first
#'   (alphabetically, for determinism) and are also reported via
#'   \code{attr(., "ambiguous")}, a named list of all candidate current
#'   symbols for each flagged old symbol.
#' @examples
#' \donttest{
#'   hgncTable <- getHgncGeneTable()
#'   symbolMap <- buildHgncSymbolMap(hgncTable)
#'   length(attr(symbolMap, "ambiguous"))
#' }
#' @seealso \code{\link{getHgncGeneTable}}, \code{\link{normalizeGeneSymbols}}
#' @export
buildHgncSymbolMap <- function(hgncTable, warn = TRUE) {
    stopifnot(is.data.frame(hgncTable),
             all(c("symbol", "prev_symbol", "alias_symbol") %in% names(hgncTable)))
    pairs <- Map(function(sym, prev, alias) {
        old <- c(prev, alias)
        old <- old[nzchar(old)]
        if (length(old) == 0L) return(NULL)
        data.frame(old = old, current = sym, stringsAsFactors = FALSE)
    }, hgncTable$symbol, hgncTable$prev_symbol, hgncTable$alias_symbol)
    pairs <- do.call(rbind, pairs)

    if (is.null(pairs)) {
        map <- character(0)
        attr(map, "ambiguous") <- list()
        return(map)
    }
    agg <- split(pairs$current, pairs$old)
    ambiguous <- agg[lengths(lapply(agg, unique)) > 1L]
    if (isTRUE(warn) && length(ambiguous)) {
        warning(length(ambiguous), " old symbol(s) map to more than one current ",
                "symbol; keeping the first (alphabetical) - see attr(., \"ambiguous\").",
                call. = FALSE)
    }
    map <- vapply(agg, function(x) sort(unique(x))[1], character(1))
    attr(map, "ambiguous") <- lapply(ambiguous, unique)
    map
}

#' Normalize gene symbols to their current HGNC-approved form
#'
#' Passes already-current symbols straight through; resolves anything
#' else via a prev/alias-symbol map (see \code{\link{buildHgncSymbolMap}});
#' anything still unresolved becomes \code{NA}, with the unmapped input
#' symbols reported via \code{attr(., "unmapped")} rather than silently
#' dropped.
#'
#' @param symbols character vector of gene symbols to normalize.
#' @param symbolMap named character vector as returned by
#'   \code{\link{buildHgncSymbolMap}}; built from \code{hgncTable} if not
#'   supplied.
#' @param hgncTable data.frame as returned by \code{\link{getHgncGeneTable}};
#'   used to build \code{symbolMap} (and to recognise already-current
#'   symbols) when \code{symbolMap} is not supplied directly. Fetched
#'   automatically if neither argument is given.
#' @return character vector, same length/order as \code{symbols}. Symbols
#'   that could not be resolved at all are \code{NA} and are listed in
#'   \code{attr(., "unmapped")}. Symbols that resolve to more than one
#'   current symbol keep the first (alphabetically) and are listed in
#'   \code{attr(., "ambiguous")}, with a \code{warning()}. Only the
#'   symbols given are considered: a symbol that is already current is
#'   passed straight through and never counts as ambiguous, so translating
#'   unaffected symbols is silent.
#' @examples
#' \donttest{
#'   hgncTable <- getHgncGeneTable()
#'   normalizeGeneSymbols(c("FGFR1", "ABL"), hgncTable = hgncTable)
#' }
#' @seealso \code{\link{buildHgncSymbolMap}}, \code{\link{getHgncGeneTable}}
#' @export
normalizeGeneSymbols <- function(symbols, symbolMap = NULL, hgncTable = NULL) {
    stopifnot(is.character(symbols))
    if (is.null(hgncTable) && is.null(symbolMap)) hgncTable <- getHgncGeneTable()
    ## warn = FALSE: how many entries the whole HGNC snapshot has that are
    ## ambiguous says nothing about the symbols being translated here. The
    ## ones that matter are reported below.
    if (is.null(symbolMap)) symbolMap <- buildHgncSymbolMap(hgncTable, warn = FALSE)

    out <- symbols
    isCurrent <- if (!is.null(hgncTable)) out %in% hgncTable$symbol else rep(FALSE, length(out))
    needsMap <- !isCurrent
    out[needsMap] <- unname(symbolMap[out[needsMap]])
    attr(out, "unmapped") <- unique(symbols[needsMap & is.na(out)])

    ## Only symbols actually being translated can hit an ambiguous entry;
    ## one already current is passed straight through untouched.
    hit <- intersect(unique(symbols[needsMap]), names(attr(symbolMap, "ambiguous")))
    attr(out, "ambiguous") <- attr(symbolMap, "ambiguous")[hit]
    if (length(hit))
        warning(length(hit), " of the symbol(s) given map to more than one ",
                "current symbol; keeping the first (alphabetical) - ",
                paste(utils::head(hit, 5), collapse = ", "),
                if (length(hit) > 5L) ", ..." else "",
                ". See attr(., \"ambiguous\").", call. = FALSE)
    out
}


## ---------------------------------------------------------------------
## Genome-wide checkpointed batch runner
## ---------------------------------------------------------------------

#' Explode an HGNC gene table's list-column into one row per (gene, id)
#' @keywords internal
#' @noRd
.hgncExplode <- function(hgncTable, col) {
    ids <- hgncTable[[col]]
    n <- lengths(ids)
    data.frame(
        hgnc_id         = rep(hgncTable$hgnc_id, n),
        symbol          = rep(hgncTable$symbol, n),
        ensembl_gene_id = rep(hgncTable$ensembl_gene_id, n),
        id              = unlist(ids, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

#' Genome-wide, checkpointed drug-target annotation build across the
#' six annotation sources, anchored on an HGNC gene table
#'
#' Loops \code{\link{getChemblDrugTarget}}/\code{\link{getDgidbDrugTarget}}/
#' \code{\link{getOpenTargetsDrugTarget}}/\code{\link{ttdTargetAnnot}}/
#' \code{\link{broadRepurposingHubAnnot}}/\code{\link{gtoPdbTargetAnnot}}
#' over every gene in \code{hgncTable} (default: all human protein-coding
#' genes from \code{\link{getHgncGeneTable}}), in checkpointed chunks
#' written to \code{outDir} as it goes - so an interrupted run resumes
#' from the last completed chunk instead of starting over. Deliberately
#' excludes ChEMBL's bioassay function and PubChem (see the file header
#' and the "Genome-Wide Master Table" vignette section for why).
#'
#' A single target-centric run already answers the reverse drug ->
#' target question too: every row already carries whichever drug
#' matched that gene, so filtering the assembled table by drug
#' identifier gives that direction for free - no separate compound-
#' centric run is needed.
#'
#' @param hgncTable data.frame from \code{\link{getHgncGeneTable}}
#'   (default: fetched automatically if not supplied).
#' @param sources character vector, any of \code{"chembl"}, \code{"dgidb"},
#'   \code{"opentargets"}, \code{"ttd"}, \code{"broad"}, \code{"gtopdb"}
#'   (default: all six).
#' @param ttdDbPath character(1) path to a local TTD SQLite (see
#'   \code{\link{buildTtdDb}}); required if \code{"ttd"} is in
#'   \code{sources}.
#' @param brhDbPath character(1) path to a local Broad Repurposing Hub
#'   SQLite (see \code{\link{buildBroadRepurposingHubDb}}); required if
#'   \code{"broad"} is in \code{sources}.
#' @param gtoPdbDbPath character(1) path to a local GtoPdb SQLite (see
#'   \code{\link{buildGtoPdbDb}}); required if \code{"gtopdb"} is in
#'   \code{sources}.
#' @param outDir character(1) directory to write checkpoint chunks and
#'   the manifest to; created if it doesn't exist.
#' @param chunkGenes integer(1) genes per checkpoint chunk (default 500)
#'   - a coarser grouping purely for resumability; each source's own
#'   function still does its own finer-grained request batching
#'   (\code{chunkSize}, passed via \code{...}) within one chunk.
#' @param rerun logical(1); if \code{FALSE} (default), chunks already
#'   recorded in \code{outDir}'s manifest are skipped - the resume
#'   behavior. Set \code{TRUE} to rebuild everything from scratch.
#' @param verbose logical(1); default \code{TRUE} (unlike this package's
#'   other \code{build*()} functions) since this is a long-running job.
#' @param ... additional arguments passed through to the ChEMBL/DGIdb/
#'   Open Targets calls (e.g. \code{fields}, \code{chunkSize}) - not
#'   forwarded to \code{ttdTargetAnnot()}/\code{broadRepurposingHubAnnot()}/
#'   \code{gtoPdbTargetAnnot()}, none of which take extra arguments.
#' @return A named list, one data.frame per successfully-built source,
#'   each row tagged with \code{hgnc_id}/\code{symbol}/\code{ensembl_gene_id}
#'   alongside that source's own native columns (not run through
#'   \code{\link{combineDrugTargets}}'s harmonization, to keep maximum
#'   information - apply that separately if wanted). Also cached via
#'   \code{\link{genConfig}}'s \code{BiocFileCache} under a name that
#'   encodes the HGNC source snapshot and build date, so a completed
#'   build is reusable without recomputation; the cache path is recorded
#'   in \code{attr(., "cachePath")}.
#' @examples
#' \donttest{
#'   ## Small illustrative run, not genome-wide - see the "Genome-Wide
#'   ## Master Table" vignette section for the full workflow (takes
#'   ## hours at full scale, not run here).
#'   hgncSmall <- getHgncGeneTable()[1:3, ]
#'   res <- buildGenomeWideDrugTargetTable(
#'     hgncTable = hgncSmall, sources = c("dgidb", "opentargets"),
#'     outDir = tempfile("dti_build_"))
#'   names(res)
#' }
#' @seealso \code{\link{getHgncGeneTable}}, \code{\link{queryDrugTargets}},
#'   \code{\link{combineDrugTargets}}
#' @export
buildGenomeWideDrugTargetTable <- function(hgncTable = NULL,
                                           sources = c("chembl", "dgidb", "opentargets",
                                                      "ttd", "broad", "gtopdb"),
                                           ttdDbPath = NULL, brhDbPath = NULL,
                                           gtoPdbDbPath = NULL, outDir, chunkGenes = 500L,
                                           rerun = FALSE, verbose = TRUE, ...) {
    if (missing(outDir) || !is.character(outDir) || length(outDir) != 1L)
        stop("'outDir' must be supplied (a single directory path).")
    sources <- match.arg(sources, c("chembl", "dgidb", "opentargets", "ttd",
                                    "broad", "gtopdb"),
                         several.ok = TRUE)
    if ("ttd" %in% sources && is.null(ttdDbPath))
        stop("'ttd' requires ttdDbPath (see buildTtdDb()).")
    if ("broad" %in% sources && is.null(brhDbPath))
        stop("'broad' requires brhDbPath (see buildBroadRepurposingHubDb()).")
    if ("gtopdb" %in% sources && is.null(gtoPdbDbPath))
        stop("'gtopdb' requires gtoPdbDbPath (see buildGtoPdbDb()).")
    if (is.null(hgncTable)) hgncTable <- getHgncGeneTable()
    stopifnot(is.data.frame(hgncTable),
             all(c("hgnc_id", "symbol", "ensembl_gene_id", "uniprot_ids") %in% names(hgncTable)))

    if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
    manifestPath <- file.path(outDir, "manifest.rds")
    manifest <- if (!rerun && file.exists(manifestPath)) readRDS(manifestPath) else character(0)

    ## Per-source key table: one row per (gene, native ID to query).
    ## ChEMBL is UniProt-accession-keyed (exploded, since a gene can have
    ## several reviewed accessions); the others are symbol-keyed (one
    ## row per gene already, no explosion needed).
    keyTables <- list()
    if ("chembl" %in% sources) keyTables$chembl <- .hgncExplode(hgncTable, "uniprot_ids")
    for (src in intersect(sources, c("dgidb", "opentargets", "ttd", "broad", "gtopdb"))) {
        keyTables[[src]] <- data.frame(
            hgnc_id = hgncTable$hgnc_id, symbol = hgncTable$symbol,
            ensembl_gene_id = hgncTable$ensembl_gene_id, id = hgncTable$symbol,
            stringsAsFactors = FALSE)
    }

    for (src in sources) {
        keyTbl <- keyTables[[src]]
        chunks <- .dtiChunk(seq_len(nrow(keyTbl)), chunkGenes)
        for (i in seq_along(chunks)) {
            chunkId <- paste0(src, "_chunk", i)
            if (chunkId %in% manifest) next
            ids <- unique(stats::na.omit(keyTbl$id[chunks[[i]]]))
            if (verbose) message("buildGenomeWideDrugTargetTable: ", src,
                                 " chunk ", i, "/", length(chunks),
                                 " (", length(ids), " ID(s))")
            result <- if (length(ids) == 0L) NULL else tryCatch({
                qb <- switch(src,
                    chembl      = list(molType = "protein", idType = "Uniprot", ids = ids),
                    dgidb       = list(molType = "gene", idType = "symbol", ids = ids),
                    opentargets = list(molType = "gene", idType = "symbol", ids = ids),
                    ttd         = list(molType = "protein", idType = "symbol", ids = ids),
                    broad       = list(molType = "protein", idType = "symbol", ids = ids),
                    gtopdb      = list(molType = "protein", idType = "symbol", ids = ids))
                switch(src,
                    chembl      = getChemblDrugTarget(qb, verbose = FALSE, ...),
                    dgidb       = getDgidbDrugTarget(qb, verbose = FALSE, ...),
                    opentargets = getOpenTargetsDrugTarget(qb, verbose = FALSE, ...),
                    ttd         = ttdTargetAnnot(qb, ttdDbPath),
                    broad       = broadRepurposingHubAnnot(qb, brhDbPath),
                    gtopdb      = gtoPdbTargetAnnot(qb, gtoPdbDbPath))
            }, error = function(e) {
                warning("buildGenomeWideDrugTargetTable: ", src, " chunk ", i,
                        " failed: ", conditionMessage(e), call. = FALSE)
                NULL
            })
            saveRDS(result, file.path(outDir, paste0(chunkId, ".rds")))
            manifest <- c(manifest, chunkId)
            saveRDS(manifest, manifestPath)
        }
    }

    ## Assemble: rbind every chunk file per source, then join gene
    ## identity back in via the (exploded) key table - QueryIDs is the
    ## accession/symbol actually queried with, which the key table
    ## already maps to hgnc_id/symbol/ensembl_gene_id.
    out <- list()
    for (src in sources) {
        chunkFiles <- list.files(outDir, pattern = paste0("^", src, "_chunk[0-9]+\\.rds$"),
                                 full.names = TRUE)
        if (length(chunkFiles) == 0L) next
        parts <- Filter(Negate(is.null), lapply(chunkFiles, readRDS))
        if (length(parts) == 0L) next
        combined <- do.call(rbind, parts)
        keyTbl <- unique(keyTables[[src]][, c("hgnc_id", "symbol", "ensembl_gene_id", "id")])
        combined <- merge(keyTbl, combined, by.x = "id", by.y = "QueryIDs", all.y = TRUE)
        names(combined)[names(combined) == "id"] <- "QueryIDs"
        rownames(combined) <- NULL
        out[[src]] <- combined
    }

    ## Derive a short tag from the HGNC source - deliberately NOT the raw
    ## hgncSource string itself. BiocFileCache::bfcquery()'s default
    ## rname matching is substring-based, not exact (confirmed live): an
    ## rname that *contains* another cache entry's exact rname as a
    ## substring will later shadow lookups for that shorter name via
    ## .getCacheFile()'s most-recent-by-create_time convention - which is
    ## exactly what happened when this used the full HGNC filename here.
    hgncSrc <- attr(hgncTable, "hgncSource") %||% "custom"
    hgncTag <- regmatches(hgncSrc, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", hgncSrc))
    if (length(hgncTag) == 0L || !nzchar(hgncTag)) hgncTag <- "custom"
    cacheName <- paste0("genomewide-drugtargets-hgnc", hgncTag, "-build",
                        format(Sys.Date(), "%Y%m%d"), ".rds")
    tmp <- tempfile(fileext = ".rds")
    saveRDS(out, tmp)
    bfc <- .getCache()
    rid <- names(bfcadd(bfc, cacheName, tmp, action = "copy"))
    file.remove(tmp)
    attr(out, "cachePath") <- bfcrpath(bfc, rids = rid)
    out
}


## ---------------------------------------------------------------------
## Assembling the build into one table
## ---------------------------------------------------------------------
## buildGenomeWideDrugTargetTable() returns one table per source, each
## with that source's own columns. Stacking them into a single table
## needs a mapping from those columns onto shared ones, which
## drugTargetColumnMap() supplies and the user can edit (columnMap.R).
##
## The rule that shapes this function: the gene identity a genome-wide
## build attaches is correct by construction - every row was queried
## *from* a known HGNC gene - and is therefore carried through untouched,
## never recomputed from what the source echoed back. Recomputing it is
## strictly worse and sometimes wrong. Measured on a real four-source
## build, re-deriving hgnc_id from each source's own identifiers recovers
## only 497 of 2000 ChEMBL rows (many UniProt accessions name more than
## one gene, so the honest answer is NA), and for TTD it produces 98 rows
## keyed to the wrong gene entirely: rows queried as ADRA1A come back as
## HGNC:280, which is ADRA1D, because ADRA1A is also a previous symbol of
## that gene and the tie-break picks alphabetically.

#' Append a genome-wide build's per-source tables into one table
#'
#' \code{\link{buildGenomeWideDrugTargetTable}} returns one table per
#' source, each keeping that source's own columns. This stacks them into
#' a single table: columns holding the same content under different names
#' are aligned onto one shared column, and columns belonging to only one
#' source are carried through as their own columns, empty for every other
#' source's rows.
#'
#' Which columns hold the same content is decided by \code{colMap}, an
#' ordinary data.frame from \code{\link{drugTargetColumnMap}} that you can
#' edit to change the alignment or to add groups of your own.
#'
#' Each row keeps the HGNC gene it was queried from, in \code{hgnc_id},
#' \code{symbol} and \code{ensembl_gene_id}. Those are taken from the
#' build unchanged. \code{symbol} is HGNC's current symbol for that gene;
#' \code{gene_symbol} is what the source itself called it, and the two
#' can differ where a source reports an older symbol.
#'
#' Some columns share a concept without sharing a vocabulary.
#' \code{max_phase} is the clearest case: ChEMBL reports a number from
#' -1 to 4, Open Targets a token such as \code{"PHASE_3"}, TTD and the
#' Broad Repurposing Hub a phrase such as \code{"Approved"}. Read that
#' column together with \code{source}.
#'
#' @param results a named list of per-source data.frames from
#'   \code{\link{buildGenomeWideDrugTargetTable}}, or any similarly-shaped
#'   list such as \code{\link{queryDrugTargets}}' output.
#' @param colMap the column mapping to apply; default
#'   \code{\link{drugTargetColumnMap}()}, the curated mapping the package
#'   ships. Rows with \code{active = FALSE} are ignored.
#' @param native logical(1); if \code{TRUE} (default), columns belonging
#'   to a single source are carried through as \code{source.column}
#'   (\code{ttd.Smiles}, \code{opentargets.disease_id}) rather than
#'   dropped, so nothing in the build is lost. Set \code{FALSE} for just
#'   the shared columns.
#' @param keys logical(1); if \code{TRUE} (default), make sure
#'   \code{hgnc_id} and \code{compound_chembl_id} are present, deriving
#'   with \code{\link{addCommonIds}} only what is missing. Values already
#'   in \code{results} are never overwritten.
#' @param hgncTable optional HGNC table from
#'   \code{\link{getHgncGeneTable}}, passed to \code{\link{addCommonIds}}
#'   when keys have to be derived.
#' @param verbose logical(1); report per-source row counts and any
#'   columns the mapping does not cover (default \code{FALSE}).
#' @return A single \code{data.frame}: the gene identity columns present
#'   in \code{results}, then \code{source}, then one column per shared
#'   group in \code{colMap}, then \code{compound_chembl_id}, then each
#'   source's own remaining columns when \code{native = TRUE}.
#' @examples
#' \donttest{
#'   res <- buildGenomeWideDrugTargetTable(
#'     hgncTable = getHgncGeneTable()[1:3, ],
#'     sources = c("dgidb", "opentargets"), outDir = tempfile("dti_"))
#'
#'   ## Shared columns only
#'   combineGenomeWideDrugTargets(res, native = FALSE)
#'
#'   ## Everything, with each source's own columns kept
#'   tbl <- combineGenomeWideDrugTargets(res)
#'   table(tbl$source)
#' }
#' @seealso \code{\link{buildGenomeWideDrugTargetTable}},
#'   \code{\link{drugTargetColumnMap}}, \code{\link{addCommonIds}},
#'   \code{\link{mergeDrugTargets}}
#' @export
combineGenomeWideDrugTargets <- function(results,
                                         colMap = drugTargetColumnMap(),
                                         native = TRUE, keys = TRUE,
                                         hgncTable = NULL, verbose = FALSE) {
    if (!is.list(results) || is.data.frame(results))
        stop("'results' must be a named list of per-source data.frames, as ",
             "returned by buildGenomeWideDrugTargetTable() or ",
             "queryDrugTargets().")
    if (length(results) == 0L) return(data.frame())
    unsupported <- setdiff(names(results), names(.dtiCombineSourceLabel))
    if (length(unsupported))
        stop("combineGenomeWideDrugTargets() does not recognise source(s): ",
             paste(unsupported, collapse = ", "), ". Expected any of: ",
             paste(names(.dtiCombineSourceLabel), collapse = ", "), ".")

    colMap <- .dtiValidateColumnMap(colMap)

    ## Derive only the keys that are actually missing. addCommonIds()
    ## leaves an hgnc_id that is already there alone, so a genome-wide
    ## build keeps the gene identity it established.
    if (isTRUE(keys)) {
        lacking <- vapply(results, function(d)
            is.data.frame(d) && nrow(d) > 0L &&
            !all(.dtiKeyOwnedCols %in% names(d)), logical(1))
        if (any(lacking))
            results <- addCommonIds(results, hgncTable = hgncTable,
                                    verbose = verbose)
    }

    nonEmpty <- vapply(results, function(d)
        is.data.frame(d) && nrow(d) > 0L, logical(1))
    results <- results[nonEmpty]
    if (length(results) == 0L) return(data.frame())

    ## Layout, left to right.
    idCols <- .dtiIdentityCols[vapply(.dtiIdentityCols, function(cl)
        any(vapply(results, function(d) cl %in% names(d), logical(1))), logical(1))]
    canonCols <- unique(colMap$canonical)
    keyCols <- if (isTRUE(keys)) "compound_chembl_id" else character(0)
    keyCols <- keyCols[vapply(keyCols, function(cl)
        any(vapply(results, function(d) cl %in% names(d), logical(1))), logical(1))]

    ## Native columns are whatever a source carries that neither the
    ## identity block, the derived keys, nor an active mapping row
    ## already accounts for.
    nativeBySrc <- lapply(names(results), function(src) {
        if (!isTRUE(native)) return(character(0))
        mapped <- colMap$column[colMap$source == src]
        setdiff(names(results[[src]]),
                c(.dtiIdentityCols, .dtiCommonIdCols, mapped))
    })
    names(nativeBySrc) <- names(results)
    nativeCols <- unlist(lapply(names(results), function(src)
        if (length(nativeBySrc[[src]])) paste(src, nativeBySrc[[src]], sep = ".")
        else NULL), use.names = FALSE)
    nativeCols <- unique(c(nativeCols))

    outCols <- c(idCols, "source", canonCols, keyCols, nativeCols)

    blocks <- lapply(names(results), function(src) {
        df <- results[[src]]
        n <- nrow(df)
        out <- data.frame(row.names = seq_len(n))
        ## Identity and derived keys: verbatim, never recomputed here.
        for (cl in c(idCols, keyCols))
            out[[cl]] <- if (cl %in% names(df)) as.character(df[[cl]]) else NA_character_
        out$source <- .dtiCombineSourceLabel[[src]]
        ## Shared columns: as.character() throughout, since sources that
        ## share a concept do not always share its type (Max_Phase is a
        ## number in ChEMBL and a phrase in TTD).
        for (cl in canonCols) {
            srcCol <- colMap$column[colMap$canonical == cl & colMap$source == src]
            out[[cl]] <- if (length(srcCol) == 1L && srcCol %in% names(df))
                as.character(df[[srcCol]]) else NA_character_
        }
        ## Source-unique columns keep their own type - each belongs to
        ## exactly one source, so nothing has to be reconciled.
        for (cl in nativeCols) out[[cl]] <- NA
        for (cl in nativeBySrc[[src]]) out[[paste(src, cl, sep = ".")]] <- df[[cl]]
        if (verbose) {
            unmapped <- setdiff(names(df), c(.dtiIdentityCols, .dtiCommonIdCols,
                                             colMap$column[colMap$source == src]))
            message("combineGenomeWideDrugTargets: ", src, " - ", n, " row(s)",
                    if (length(unmapped))
                        paste0(", ", length(unmapped), " column(s) outside the ",
                               "mapping: ", paste(unmapped, collapse = ", "),
                               if (isTRUE(native)) " (carried through)"
                               else " (dropped)")
                    else "")
        }
        out[, outCols, drop = FALSE]
    })

    out <- do.call(rbind, blocks)
    rownames(out) <- NULL
    out
}
