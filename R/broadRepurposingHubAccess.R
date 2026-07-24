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
##  Known operational caveat: as of this writing repo-hub.broadinstitute.org's
##  TLS certificate chain is served without its InCommon intermediate
##  certificate, which can cause certificate verification failures on
##  machines/OSes whose trust store hasn't independently cached that
##  intermediate (most browsers chase the AIA "CA Issuers" URL
##  automatically; curl/libcurl/R's default download methods generally
##  do not). This is a server-side misconfiguration, not something to
##  paper over by disabling verification - so downloadBroadRepurposingHub()
##  instead retries once, on any download failure, with a CA bundle
##  extended to include *only* that one legitimate, publicly-issued
##  intermediate (see .brhWithSupplementalCa() below) - full chain
##  verification still happens, it just has the one certificate Broad's
##  own server should have sent but doesn't.
## =====================================================================


## ---------------------------------------------------------------------
## Internal infrastructure
## ---------------------------------------------------------------------

#' Endpoint registry for the Broad Repurposing Hub flat files
#' @keywords internal
#' @noRd
.brhEndpoints <- function() {
    list(
        drug   = "https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20200324.txt",
        sample = "https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20240610.txt"
    )
}

## InCommon RSA OV SSL CA 3 - the intermediate repo-hub.broadinstitute.org
## fails to send during its TLS handshake (see file header). Retrieved
## once from the cert's own "CA Issuers" AIA URI
## (http://crt.sectigo.com/InCommonRSAOVSSLCA3.crt) and embedded here so
## no extra live network call is needed just to work around this.
.brhIncommonIntermediatePem <- "-----BEGIN CERTIFICATE-----
MIIGIzCCBAugAwIBAgIRAJa22zsNXLm6Xd6KrItt/ncwDQYJKoZIhvcNAQEMBQAw
XzELMAkGA1UEBhMCR0IxGDAWBgNVBAoTD1NlY3RpZ28gTGltaXRlZDE2MDQGA1UE
AxMtU2VjdGlnbyBQdWJsaWMgU2VydmVyIEF1dGhlbnRpY2F0aW9uIFJvb3QgUjQ2
MB4XDTI1MTEwNjAwMDAwMFoXDTM1MTEwNTIzNTk1OVowSDELMAkGA1UEBhMCVVMx
FjAUBgNVBAoTDUluQ29tbW9uLCBMTEMxITAfBgNVBAMTGEluQ29tbW9uIFJTQSBP
ViBTU0wgQ0EgMzCCAaIwDQYJKoZIhvcNAQEBBQADggGPADCCAYoCggGBAInzD7j/
Ja1OZOvyIIe2hFOdDrois8Iiuyh+RtSKaKyQAvSRdG1b0Iz+fxOZaNPlM2RCTa9N
Ar/bs9Tts4RXTDCuLJfCPPwbRtSZMvBrZVpcPU3xVBbTUHTsYZ+SmlzB+qIwEJV6
TU8vEsdqosCwA/iOXewiRmUf5FxU2WoU4nD8iVFhu/p6h6YmI+AgswZ4lwZdNKW5
9cTvpuY8VefEWHuwvSQlzekLLBqiFJhlCu8dNrBsahT07sjMHVZVHU8Biss3bX04
FTzkDzv5eZ/U2LFA0rV2QzLpeLtIsMsXhlrEmuT4g6cbJJ3ZfWGHX77jnCIczshi
taD1BiTA9PRv9JWW6xQ+cGfHRMyHWTBNhdQI22N9UO65R+6ddwEWupViEGRUuO/O
ZbTtBSkoEBejHBfI3/BnnhLpKXTGC5N20om7nQ2UqOcgOpewiE9P+DnMrnUsqp9e
IS6NgkCCCuhS6eoHS3DwJouIK1T5CE/4xSrEuU0QTFRJfyqlIaMbKMe+awIDAQAB
o4IBbzCCAWswHwYDVR0jBBgwFoAUVnNYZJX5khqwEioEYnmhQBWIIUkwHQYDVR0O
BBYEFNoiNz/l03Ta2Xk+0XJt1ZNLIDevMA4GA1UdDwEB/wQEAwIBhjASBgNVHRMB
Af8ECDAGAQH/AgEAMBMGA1UdJQQMMAoGCCsGAQUFBwMBMBMGA1UdIAQMMAowCAYG
Z4EMAQICMFQGA1UdHwRNMEswSaBHoEWGQ2h0dHA6Ly9jcmwuc2VjdGlnby5jb20v
U2VjdGlnb1B1YmxpY1NlcnZlckF1dGhlbnRpY2F0aW9uUm9vdFI0Ni5jcmwwgYQG
CCsGAQUFBwEBBHgwdjBPBggrBgEFBQcwAoZDaHR0cDovL2NydC5zZWN0aWdvLmNv
bS9TZWN0aWdvUHVibGljU2VydmVyQXV0aGVudGljYXRpb25Sb290UjQ2LnA3YzAj
BggrBgEFBQcwAYYXaHR0cDovL29jc3Auc2VjdGlnby5jb20wDQYJKoZIhvcNAQEM
BQADggIBADoIPZD+zZzMmsZaUIc4WiV5NwHbB5nnmvSaDas20GSsyRiSnQVUwCT6
RzJJhPGJnoIHL7uYyjZDnYrB4MOL/1c0g+7BFmDY+0/csUwdHlouTOrj17T3nyrR
JjEg3/bY16ojl91ji4g4XbvB7L2tKkK2kF+dcIECRbcw+vE2gLSv7wcl78+m0jjb
3nw8Z+bs/R/W7C8kn+6bfgRrI2NGhd2wnJ579xMLUoodj/L2sXokw0/jiDrWgGAd
MijIVvVFSTaI08/8LReuluxbFCvftNBwiBVZm7UMV3hwZ97dqW+Tq4+Lh9GrbnO/
tMMSVSib+KKRMDYh5HGfjmmW9UTmaTL23oP7XPNuuOwqy0Z4RGqMWd8JYvNe0XFw
7rs73SJJN3zccKZvznDSaCzJBCjT3i3JrbbSW1cZKn0RFHhkFZmtP63HiXAo+G0n
Z7C/INTWdQcy9tJ/tYfC/ZVjap7R2C8s7XO2PR3oBvATgaRFIPF7q7+kyN0qmwfl
7Kt0O+dFw88fvjKfMIMQtnNqY7bWpGWg0XstM1L4kwi1FFNElCjFJR3rN5kpf5n1
7LqhsLc296cnFIwf/Yu9FaPkOJiGga6FewBTGqcWj8lL2KgHn95hdky0FzgIKRYw
loK8Ee0Y9wee43mVaHEdRD11fUQZYXsnXIFFtVXsIFH43LTiYfgN
-----END CERTIFICATE-----
"

## Candidate default CA bundle file locations across common platforms -
## used only as a *base* to append the missing intermediate to, never
## as a replacement for the OS's own trust decisions. The first one
## that exists is used; if none exist (unrecognized platform, or a
## macOS build that verifies via Keychain rather than a PEM file),
## .brhWithSupplementalCa() runs its expression unmodified rather than
## guessing further, so the plain code path's exact behavior is the
## only possible fallback - never a weaker one.
.brhSystemCaBundleCandidates <- function() {
    c(
        Sys.getenv("CURL_CA_BUNDLE", NA_character_),
        Sys.getenv("SSL_CERT_FILE", NA_character_),
        "/etc/ssl/certs/ca-certificates.crt",   # Debian/Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",     # RHEL/Fedora/CentOS
        "/etc/ssl/cert.pem",                    # Alpine; some macOS
        "/usr/local/etc/openssl/cert.pem",      # Homebrew OpenSSL, Intel Mac
        "/opt/homebrew/etc/openssl@3/cert.pem"  # Homebrew OpenSSL, Apple Silicon
    )
}

#' Retry an expression with a supplemental CA bundle
#'
#' Runs \code{expr} with \code{CURL_CA_BUNDLE} pointed at [system default
#' bundle + the missing InCommon intermediate], restoring the prior
#' value afterward. Falls back to running \code{expr} completely
#' unmodified if no base bundle can be located - this only ever *adds*
#' one legitimate, publicly-issued certificate on top of whatever the
#' platform already trusts; it never disables or weakens verification.
#' @keywords internal
#' @noRd
.brhWithSupplementalCa <- function(expr) {
    candidates <- .brhSystemCaBundleCandidates()
    hits <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
    if (length(hits) == 0L) return(expr)
    base <- hits[1]

    tmp <- tempfile(fileext = ".pem")
    writeLines(c(readLines(base, warn = FALSE), .brhIncommonIntermediatePem), tmp)

    old <- Sys.getenv("CURL_CA_BUNDLE", unset = NA_character_)
    Sys.setenv(CURL_CA_BUNDLE = tmp)
    on.exit({
        if (is.na(old)) Sys.unsetenv("CURL_CA_BUNDLE") else Sys.setenv(CURL_CA_BUNDLE = old)
        unlink(tmp)
    })
    expr
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
#' @noRd
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
#' \code{repo-hub.broadinstitute.org} has been observed to serve an
#' incomplete TLS certificate chain (see this file's header comment) -
#' if the first download attempt fails for any reason, this function
#' retries once with a CA bundle extended to include the one
#' legitimate, publicly-issued intermediate certificate the server
#' itself omits (\code{.brhWithSupplementalCa()}); a genuine failure
#' (network down, file moved, etc.) still surfaces as an error after
#' that retry.
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
            paths[[nm]] <- tryCatch(
                suppressWarnings(.downloadFile(ep[[nm]], fname)),
                error = function(e) .brhWithSupplementalCa(.downloadFile(ep[[nm]], fname))
            )
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
#' @noRd
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


#' Derive a compound-level (one row per \code{pert_iname}) structure
#' lookup from the sample table: \code{smiles}/\code{InChIKey}/
#' \code{pubchem_cid}.
#'
#' \code{sample} is one row per physical sample (\code{broad_id}), not
#' one row per compound - a compound with several samples/lots
#' otherwise fans out any merge keyed on \code{pert_iname} (confirmed
#' live, 2026-07-22: this was the root cause of a real bug, see
#' \code{\link{buildBroadRepurposingHubDb}}). The overwhelming majority
#' of compounds (verified live: ~98%) have exactly one structure across
#' all their samples, collapsing to a single row with no ambiguity.
#' Before comparing structures across samples, incomplete records (any
#' of the three columns \code{NA}) are dropped whenever at least one
#' complete record exists for that compound - a sample missing e.g.
#' \code{InChIKey} while another sample of the same compound has it
#' fully populated is missing metadata, not a second structure
#' (confirmed live: naively comparing raw distinct combinations flagged
#' 312 compounds as multi-structure, but 182 of those were exactly this
#' NA-driven false positive - only ~130-132 are genuinely
#' multi-structure once complete records are preferred). Of those, a
#' small minority (~2% of all compounds, as of the 2025-08-18 release)
#' have a \code{pert_iname} display name that actually covers more than
#' one distinct structure - different salt forms or stereoisomers
#' sharing one name, a real property of the source data, not a parsing
#' artifact. For those, every distinct structure is packed into one
#' \code{"; "}-joined string per column (same convention as
#' \code{ttdAccess.R}'s \code{Indication} packing), with the i-th
#' segment aligned across all three columns so a value from one
#' structure is never paired with another's, and
#' \code{structure_ambiguous} is set \code{TRUE} - callers doing
#' structure-sensitive work should treat \code{pert_iname} as a display
#' label hiding real chemical multiplicity for those rows, not a
#' checkable molecule identity; the honest key there is
#' \code{InChIKey}, not name.
#' @keywords internal
#' @noRd
.brhCompoundStructure <- function(sample) {
    structCols <- c("smiles", "InChIKey", "pubchem_cid")
    distinct <- unique(sample[, c("pert_iname", structCols)])
    grouped <- split(distinct[structCols], distinct$pert_iname)
    out <- do.call(rbind, lapply(names(grouped), function(nm) {
        d <- grouped[[nm]]
        if (nrow(d) > 1L) {
            ## Prefer complete-case rows over incomplete ones first -
            ## a sample missing e.g. InChIKey while another sample of
            ## the same compound has it fully populated is missing
            ## metadata, not a second structure (confirmed live:
            ## 182/312 of the naive multi-combo count was exactly this
            ## NA-driven false positive; only 130-132 are genuinely
            ## multi-structure once complete rows are preferred).
            complete <- d[stats::complete.cases(d), , drop = FALSE]
            if (nrow(complete) >= 1L) d <- unique(complete)
        }
        if (nrow(d) == 1L) {
            data.frame(pert_iname = nm, smiles = d$smiles, InChIKey = d$InChIKey,
                      pubchem_cid = d$pubchem_cid, structure_ambiguous = FALSE,
                      stringsAsFactors = FALSE)
        } else {
            data.frame(pert_iname = nm,
                      smiles      = paste(d$smiles, collapse = "; "),
                      InChIKey    = paste(d$InChIKey, collapse = "; "),
                      pubchem_cid = paste(d$pubchem_cid, collapse = "; "),
                      structure_ambiguous = TRUE, stringsAsFactors = FALSE)
        }
    }))
    rownames(out) <- NULL
    out
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
#' \code{\link{downloadBroadRepurposingHub}}) and writes **two** tables
#' to a local SQLite file, each kept at its own true grain rather than
#' flattened into one (see "Two-table design" below):
#' \code{broad_interactions} (one row per \code{(pert_iname,
#' target_gene)} drug-target edge) and \code{broad_samples} (physical-
#' sample/QC metadata, one or more rows per \code{broad_id}). The file
#' is cached via BiocFileCache under a name that includes the source
#' files' own \code{"!File_date"} (e.g. \code{broad_repurposing_20250818.db}),
#' so rebuilding is a no-op until the Hub actually republishes. As with
#' \code{\link{buildTtdDb}}, the package itself never ships or
#' redistributes this file - it is built into the caller's own local
#' cache the first time this is run.
#'
#' @section Two-table design (2026-07-22):
#' Originally a single flat \code{broad_interactions} table, mirroring
#' \code{\link{buildTtdDb}}'s one-table shape - but unlike TTD, the
#' Repurposing Hub's sample-level file (\code{sample}, one row per
#' physical vial/lot, \code{broad_id}) genuinely has a *different*
#' grain than the drug-target edge the interaction table is meant to
#' represent: a single compound (\code{pert_iname}) routinely has
#' several physical samples (verified live: 67% of compounds have more
#' than one \code{broad_id}, up to 8), and each sample can itself carry
#' several repeat QC/purity readings. A one-table design merging
#' \code{sample} straight into the exploded target table (via
#' \code{pert_iname}, deduplicated only by
#' \code{sample[!duplicated(sample), ]}, which drops exact duplicate
#' *rows*, not duplicate *join keys*) confirmed to multiply rows in
#' practice - a real, shipped bug: 3.57x overall row inflation
#' (17,996 true drug-target edges vs. 64,336 rows), with one
#' compound-target pair duplicated 175x, driven by repeat purity
#' readings for a single \code{broad_id} that survived the row-level
#' dedup untouched. Fixed by never merging \code{sample} into the
#' interaction table at all: \code{broad_interactions} draws its
#' compound-level structure fields from
#' \code{.brhCompoundStructure} instead (a proper one-row-per-
#' compound reduction of \code{sample}, packing the rare cases where a
#' name genuinely covers multiple structures rather than picking one
#' arbitrarily), and everything else \code{sample} carries (purity,
#' vendor, catalog_no, \code{qc_incompatible}, etc.) - which was never
#' actually a property of a drug-target edge to begin with, only of a
#' physical vial - moves to \code{broad_samples}, kept at its own
#' natural grain. \code{broad_interactions}' grain is asserted, not
#' just assumed, via \code{.assertUniqueKey} immediately after
#' the reduction, so a future Repurposing Hub release that breaks this
#' assumption fails the build loudly instead of silently reinflating
#' results. \code{idType = "broad_id"} queries
#' (\code{\link{broadRepurposingHubAnnot}}) now resolve \code{broad_id
#' -> pert_iname} via \code{broad_samples} first, then dispatch into
#' \code{broad_interactions} - which, given the 67% multi-sample rate,
#' is also a correctness improvement for that query path, not just a
#' side effect of the split.
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

    ## broad_interactions: one row per (pert_iname, target_gene) drug-
    ## target edge - see "Two-table design" above for why structure
    ## comes from a compound-level reduction of `sample`
    ## (.brhCompoundStructure()) rather than a raw merge against it.
    structure <- .brhCompoundStructure(sample)
    exploded  <- .brhExplodeTargets(drug)
    interactions <- merge(exploded, structure, by = "pert_iname", all.x = TRUE)
    interactions <- interactions[, c("target_gene", "pert_iname", "clinical_phase",
                                     "moa", "disease_area", "indication",
                                     "smiles", "InChIKey", "pubchem_cid",
                                     "structure_ambiguous")]
    .assertUniqueKey(interactions, c("pert_iname", "target_gene"), "broad_interactions")

    ## broad_samples: physical-sample/QC-level metadata (purity, vendor,
    ## etc. - never actually a property of a drug-target edge), kept at
    ## its own natural grain rather than forced to one row per broad_id
    ## - a broad_id can legitimately have several rows here (e.g.
    ## repeat purity checks on the same vial), which is faithful to the
    ## source data, not a bug, since this table is never merged back
    ## into broad_interactions. Used to resolve idType="broad_id"
    ## queries (broadRepurposingHubAnnot()) and as a QC reference table
    ## in its own right.
    samples <- sample[, c("broad_id", "pert_iname", "qc_incompatible", "purity",
                          "vendor", "catalog_no", "vendor_name", "expected_mass",
                          "deprecated_broad_id")]
    samples <- samples[!duplicated(samples), ]

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "broad_interactions", interactions, overwrite = TRUE)
    dbWriteTable(con, "broad_samples", samples, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_brh_target    ON broad_interactions (target_gene)")
    dbExecute(con, "CREATE INDEX idx_brh_pert      ON broad_interactions (pert_iname)")
    dbExecute(con, "CREATE INDEX idx_brhs_broad_id ON broad_samples (broad_id)")
    dbExecute(con, "CREATE INDEX idx_brhs_pert     ON broad_samples (pert_iname)")
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
#' \code{target_gene = NA} (see \code{.brhExplodeTargets}), so it
#' remains reachable by name/broad_id even though it can never surface
#' via a gene-side query.
#'
#' \code{idType = "broad_id"} does not query \code{broad_interactions}
#' directly - since the 2026-07-22 two-table split (see
#' \code{\link{buildBroadRepurposingHubDb}}), \code{broad_id} lives only
#' in \code{broad_samples}. Each queried \code{broad_id} is first
#' resolved to its compound name(s) there, then dispatched exactly like
#' an \code{idType = "name"} lookup, with \code{QueryIDs} kept tagged to
#' the original \code{broad_id} rather than the resolved name. A
#' \code{broad_id} resolves to exactly one compound in the overwhelming
#' majority of cases; a small number (10 as of the 2025-08-18 release)
#' resolve to more than one due to genuine naming inconsistencies in
#' Broad's own source file (e.g. \code{"prednisolone acetate"} vs.
#' \code{"prednisolone-acetate"} sharing one \code{broad_id}) - handled
#' as a real one-to-many resolution (all matching compounds' full
#' target sets are returned) rather than collapsed to one or errored on.
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
#'   \code{moa}, \code{disease_area}, \code{indication}, \code{smiles},
#'   \code{InChIKey}, \code{pubchem_cid}, \code{structure_ambiguous}
#'   (with the default \code{fields = "core"}), or a subset when
#'   \code{fields} requests specific columns. \code{structure_ambiguous
#'   = TRUE} means \code{pert_iname} covers more than one distinct
#'   structure for that row (packed \code{"; "}-joined into
#'   \code{smiles}/\code{InChIKey}/\code{pubchem_cid}, aligned segment-
#'   by-segment across the three columns) - a real name collision in
#'   the source data (different salts/stereoisomers sharing a display
#'   name), not a parsing artifact; structure-sensitive work should key
#'   on \code{InChIKey}, not \code{pert_iname}, for those rows. Sample-
#'   level fields (\code{purity}, \code{vendor}, \code{catalog_no},
#'   \code{qc_incompatible}, \code{expected_mass},
#'   \code{deprecated_broad_id}) moved to \code{broad_samples} in the
#'   same SQLite file - query it directly (keyed on \code{broad_id})
#'   for that data; it is not returned here.
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
               broad_id = "pert_iname",
               stop("idType for molType='cmp' must be one of: ",
                    "'name', 'broad_id'"))
    } else {
        stop("molType must be 'protein'/'gene' or 'cmp'")
    }

    con <- dbConnect(SQLite(), brhDbPath)
    on.exit(dbDisconnect(con))

    if (queryBy$idType == "broad_id") {
        ## broad_id lives only in broad_samples since the two-table
        ## split - resolve to compound name(s) first (exact-match, not
        ## case-folded), then dispatch like a "name" lookup below, but
        ## build index_list from the ORIGINAL broad_id tokens via this
        ## resolution map rather than direct token comparison.
        rIdvec <- paste0("('", paste(gsub("'", "''", queryBy$ids, fixed = TRUE), collapse = "', '"), "')")
        resolveMap <- dbGetQuery(con, paste0(
            "SELECT broad_id, pert_iname FROM broad_samples WHERE broad_id IN ", rIdvec))
        lookupNames <- unique(resolveMap$pert_iname)
    } else {
        lookupNames <- queryBy$ids
    }

    caseInsensitive <- TRUE
    ids <- toupper(lookupNames)
    resultDF <- if (length(ids) == 0L) {
        ## no broad_id resolved to anything - still need a correctly-
        ## typed 0-row frame so downstream NA-padding works the same
        ## as any other all-miss query.
        dbGetQuery(con, "SELECT * FROM broad_interactions LIMIT 0")
    } else {
        idvec <- paste0("('", paste(gsub("'", "''", ids, fixed = TRUE), collapse = "', '"), "')")
        colExpr <- paste0("UPPER(", col, ")")
        dbGetQuery(con, paste0("SELECT * FROM broad_interactions WHERE ", colExpr, " IN ", idvec))
    }

    ## Tag every row with the original query token it matched, and NA-pad
    ## any query ID that matched nothing - mirrors ttdTargetAnnot()'s
    ## QueryIDs convention exactly (same Inf-index trick: an unmatched ID
    ## gets rowid Inf, and indexing a data.frame with Inf yields a row of
    ## NAs).
    cmpFun <- toupper
    index_list <- if (queryBy$idType == "broad_id") {
        lapply(queryBy$ids, function(origId) {
            resolvedNames <- resolveMap$pert_iname[resolveMap$broad_id == origId]
            if (length(resolvedNames) == 0L) return(integer(0))
            which(cmpFun(resultDF[[col]]) %in% cmpFun(resolvedNames))
        })
    } else {
        lapply(queryBy$ids, function(x) which(cmpFun(resultDF[[col]]) %in% cmpFun(x)))
    }
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
#' @noRd
.dtiBroadAllCols <- c("QueryIDs", "target_gene", "pert_iname", "clinical_phase",
                      "moa", "disease_area", "indication", "smiles", "InChIKey",
                      "pubchem_cid", "structure_ambiguous")
