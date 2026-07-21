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
##    P1-03-TTD_crossmatching.txt     drug cross-reference IDs (PubChem
##                                     CID/SID, CAS, ChEBI) - confirmed by
##                                     reading the file body that this
##                                     release carries no ChEMBL_ID/
##                                     DrugBank_ID/InChIKey; those need a
##                                     separate UniChem hop, not this file
##    P1-05-Drug_disease.txt          per-drug disease indications with
##                                     ICD-11 codes and a *per-disease*
##                                     clinical status - richer than the
##                                     single Highest_status already in
##                                     ttd_interactions
##
##  TTD's own dedicated UniProt file (P2-01-TTD_uniprot_all.txt) was
##  checked and rejected as a source for a UniProt *accession*: its
##  UNIPROID field is the same mnemonic entry name already parsed here
##  (e.g. "FGFR1_HUMAN"), not an accession, and 613/4299 (14%) of its
##  records are literally "NOUNIPROTAC" (no TTD-supplied ID at all) -
##  overwhelmingly family/pathway/element-level "targets" with no single
##  accession by construction (e.g. "Fibroblast growth factor receptor
##  (FGFR)" vs the single-protein "FGFR1"), not a TTD data-quality bug.
##  Accession resolution instead goes through the package's own
##  GeneName -> UniProt REST scaffold (.resolveGeneIds(), idTranslation.R)
##  - see .ttdResolveUniprotAcc().
##
##  Every TTD flat file embeds its own release version on line 3 of its
##  header, e.g. "Version 10.1.01 (2024.01.10)" - this is TTD's real
##  upstream version (their equivalent of a ChEMBL release number) and is
##  used to name/cache the derived SQLite file, NOT an arbitrary download
##  timestamp. Confirmed live: TTD does not bump this version string
##  consistently across files - P1-05's copy stamps a different date
##  (2024.03.30) than P1-01/P1-03/P2-01's (2024.01.10) despite an
##  identical "10.1.01" version number - so the per-file (version, date)
##  is recorded separately for every downloaded file in the built
##  database's ttd_release table (see buildTtdDb()), not collapsed to one
##  build-wide date.
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
        base        = "https://ttd.idrblab.cn/files/download",
        targets     = "P1-01-TTD_target_download.txt",
        drugs       = "P1-02-TTD_drug_download.txt",
        mapping     = "P1-07-Drug-TargetMapping.xlsx",
        crossmatch  = "P1-03-TTD_crossmatching.txt",
        drugDisease = "P1-05-Drug_disease.txt"
    )
}

#' Extract TTD's own release version/date from a downloaded flat file header
#'
#' TTD stamps every flat file with "Version X.Y.Z (YYYY.MM.DD)" on line 3
#' of its header block - except the xlsx mapping file, which has no such
#' inline text header readable this way, so it always returns NA/NA.
#' Confirmed live: different files under the same nominal version number
#' can carry different dates (see the file header note above), so this
#' is deliberately per-file, not a single build-wide date.
#'
#' @param path character(1) path to a downloaded TTD file.
#' @return list with \code{version} and \code{date} character(1) each
#'   (e.g. \code{"10.1.01"}/\code{"2024.01.10"}), or \code{NA} for either
#'   if not found.
#' @keywords internal
.ttdRelease <- function(path) {
    if (grepl("\\.xlsx$", path, ignore.case = TRUE))
        return(list(version = NA_character_, date = NA_character_))
    hdr <- tryCatch(readLines(path, n = 5L, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character(0))
    hit <- grep("^Version ", hdr, value = TRUE)
    if (length(hit) == 0L) return(list(version = NA_character_, date = NA_character_))
    line <- hit[1]
    list(
        version = sub("^Version ([0-9.]+).*$", "\\1", line),
        date    = sub("^Version [0-9.]+ \\(([0-9.]+)\\).*$", "\\1", line)
    )
}

#' @rdname dot-ttdRelease
#' @keywords internal
.ttdVersion <- function(path) .ttdRelease(path)$version


## ---------------------------------------------------------------------
## Download raw TTD flat files (cached via the package's existing
## BiocFileCache helpers)
## ---------------------------------------------------------------------

#' Download the TTD flat files needed for target-drug mapping
#'
#' Downloads (or reuses previously cached copies of) the TTD target,
#' drug, drug-target-mapping, drug cross-matching and drug-disease
#' indication flat files via the package's existing
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
#'   \code{mapping}, \code{crossmatch}, \code{drugDisease}.
#' @examples
#' \donttest{
#'   paths <- downloadTTD()
#'   paths
#' }
#' @seealso \code{\link{buildTtdDb}}
#' @export
downloadTTD <- function(rerun = TRUE, config = genConfig()) {
    ep    <- .ttdEndpoints()
    files <- ep[c("targets", "drugs", "mapping", "crossmatch", "drugDisease")]
    paths <- vector("list", length(files))
    names(paths) <- names(files)
    for (nm in names(files)) {
        fname <- files[[nm]]
        if (rerun) {
            paths[[nm]] <- .downloadFile(paste0(ep$base, "/", fname), fname)
        } else {
            paths[[nm]] <- tryCatch(
                .getCacheFile(fname),
                error = function(e) {
                    stop("TTD flat file '", fname, "' is not cached yet. ",
                         "Run downloadTTD(rerun = TRUE) or buildTtdDb(rerun = TRUE) ",
                         "once to download it; subsequent rerun = FALSE calls will ",
                         "then reuse the cached copy.", call. = FALSE)
                }
            )
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

#' Parse P1-02's DRUGTYPE field into a per-DrugID lookup
#'
#' TTD does provide an explicit molecule-type field after all
#' (\code{DRUGTYPE}) - confirmed by reading the file body, not assumed
#' absent. It is not, however, a clean controlled vocabulary: ~80
#' distinct raw strings (\code{"Antibody"}, \code{"Antibody "} with a
#' trailing space, \code{"Monoclonal antibody"}, \code{"Monoclonal
#' Antibody"}, ...), and only ~70\% of drugs have it at all (29901/42939
#' in the pinned release). \code{\link{.ttdMoleculeType}} uses this as
#' its primary signal, normalized down to the 4 requested buckets, and
#' falls back to the SMILES-presence/name heuristic only for the ~30\%
#' of drugs with no \code{DRUGTYPE} value.
#' @keywords internal
.ttdParseDrugType <- function(path) {
    long <- .ttdReadLong(path, "D")
    dt <- .ttdFieldLookup(long, "DRUGTYPE")
    colnames(dt) <- c("DrugID", "DrugType")
    dt
}

#' Parse P1-03's cross-matching fields into a per-DrugID lookup
#'
#' Only \code{PUBCHCID}/\code{PUBCHSID}/\code{CASNUMBE}/\code{CHEBI_ID}
#' are extracted - confirmed by reading the file body (not just its
#' header/abbreviation block) that this release's crossmatching file
#' carries no \code{ChEMBL_ID}/\code{DrugBank_ID}/\code{InChIKey} at all;
#' getting those needs a separate hop through the package's existing
#' UniChem infrastructure (\code{unichemAccess.R}) via \code{PubChem_CID},
#' not TTD directly - deliberately left as a follow-on, not done here.
#'
#' \code{CASNUMBE}'s raw value is stored by TTD as a string like
#' \code{"CAS 68-96-2"} (confirmed not universal: 22/12163 rows have no
#' \code{"CAS "} prefix at all); the prefix is stripped here when
#' present, but the result is always kept as a character string, never
#' coerced to numeric - CAS numbers have leading zeros, hyphens, and a
#' checksum digit that make them non-arithmetic.
#'
#' \code{PUBCHCID}/\code{PUBCHSID} can each carry more than one
#' semicolon-separated candidate ID for a single drug (TTD's own format,
#' left as-is here, not reformatted). \code{PUBCHSID} (PubChem
#' *Substance* ID, not Compound) is comparatively low-value and kept only
#' because it costs nothing extra from the same file - nothing in this
#' package is built on it.
#' @keywords internal
.ttdParseDrugXref <- function(path) {
    long  <- .ttdReadLong(path, "D")
    cid   <- .ttdFieldLookup(long, "PUBCHCID")
    sid   <- .ttdFieldLookup(long, "PUBCHSID")
    cas   <- .ttdFieldLookup(long, "CASNUMBE")
    chebi <- .ttdFieldLookup(long, "CHEBI_ID")
    out <- data.frame(DrugID = unique(long$id), stringsAsFactors = FALSE)
    out$PubChem_CID <- cid$value[match(out$DrugID, cid$id)]
    out$PubChem_SID <- sid$value[match(out$DrugID, sid$id)]
    out$CAS         <- sub("^CAS\\s+", "", cas$value[match(out$DrugID, cas$id)])
    out$ChEBI_ID    <- chebi$value[match(out$DrugID, chebi$id)]
    out
}

#' Parse P1-05's per-drug INDICATI lines into one row per drug-indication
#'
#' Unlike P1-01/P1-02/P1-03 (a flat "ID <TAB> FIELD <TAB> VALUE" dump
#' handled by \code{.ttdReadLong()}), P1-05 identifies a record's drug
#' only once, via a standalone \code{"TTDDRUID <TAB> <id>"} line, with
#' every subsequent \code{INDICATI} line (0 or more) belonging to that
#' drug until the next \code{TTDDRUID} line - a genuinely stateful
#' format, hence a dedicated parser. \code{INDICATI} lines that precede
#' any \code{TTDDRUID} line are skipped by construction, since no drug ID
#' is active yet - but the file's own abbreviation-index legend block
#' itself contains a line \code{"TTDDRUID <TAB> TTD Drug ID"} that is
#' structurally indistinguishable from a real record marker (confirmed
#' live: without an explicit ID-format check this sets a bogus
#' \code{curId} of literally \code{"TTD Drug ID"} and attaches the
#' legend's own \code{INDICATI} line to it), so a \code{TTDDRUID} line is
#' only accepted as a real marker when its value matches TTD's actual
#' \code{"D" + alphanumeric}, no-space ID shape.
#'
#' The \code{ICD-11} sub-field is stored by TTD as \code{"ICD-11: <code>"}
#' (or the literal sentinel \code{"ICD-11: N.A."}/\code{"#N/A"} when
#' unavailable, confirmed live: 170/2 of 30315 rows respectively); the
#' \code{"ICD-11: "} prefix is stripped and the sentinels normalized to
#' \code{NA}.
#' @keywords internal
.ttdParseDrugIndications <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    parts <- strsplit(lines, "\t", fixed = TRUE)

    curId  <- NA_character_
    drugId <- character(0); disease <- character(0)
    icd11  <- character(0); status  <- character(0)
    for (p in parts) {
        if (length(p) < 2L) next
        if (identical(p[1], "TTDDRUID") && grepl("^D[0-9A-Za-z]+$", p[2])) {
            curId <- p[2]
        } else if (identical(p[1], "INDICATI") && length(p) >= 4L && !is.na(curId)) {
            drugId  <- c(drugId, curId)
            disease <- c(disease, p[2])
            icd11   <- c(icd11, p[3])
            status  <- c(status, p[4])
        }
    }
    icd11 <- sub("^ICD-11:\\s*", "", icd11)
    icd11[icd11 %in% c("N.A.", "#N/A", "")] <- NA_character_
    data.frame(DrugID = drugId, Disease = disease, ICD11 = icd11,
              ClinicalStatus = status, stringsAsFactors = FALSE)
}

#' Collapse \code{.ttdParseDrugIndications()}'s long table to one row per DrugID
#'
#' A drug can carry several disease indications, each with its own
#' clinical status - the whole reason to add this data is to recover a
#' per-disease status, unlike the single \code{Highest_status} already in
#' \code{ttd_interactions}. Packed into one \code{"; "}-delimited
#' \code{Indication} string per drug - mirroring the same
#' collapse-multiple-rows-into-one-string convention already used by the
#' package's ChEMBL-REST side for \code{Mesh_Indication} (see
#' \code{.meshFor()} in \code{apiAccess.R}) - rather than exploding
#' \code{ttd_interactions}' one-row-per-edge grain, which would multiply
#' every edge by its drug's indication count for no benefit to callers
#' who only want \code{Highest_status}. Each packed entry is
#' \code{"<disease> [<ICD-11>]: <clinical status>"} (or without the
#' bracketed code when ICD-11 is unavailable) - disease, code and status
#' are kept together per entry rather than split across three
#' separately-packed columns, so a caller never has to assume positional
#' alignment across columns to know which status goes with which disease.
#' @keywords internal
.ttdPackIndications <- function(ind) {
    if (nrow(ind) == 0L)
        return(data.frame(DrugID = character(0), Indication = character(0),
                          stringsAsFactors = FALSE))
    entry <- ifelse(is.na(ind$ICD11),
                    paste0(ind$Disease, ": ", ind$ClinicalStatus),
                    paste0(ind$Disease, " [", ind$ICD11, "]: ", ind$ClinicalStatus))
    packed <- tapply(entry, ind$DrugID, paste, collapse = "; ")
    data.frame(DrugID = names(packed), Indication = unname(packed),
              stringsAsFactors = FALSE)
}

#' Resolve TTD's mnemonic GeneName to a UniProt accession
#'
#' TTD's own \code{Uniprot} field (from \code{.ttdParseTargets()}) is a
#' mnemonic entry name (e.g. \code{"FGFR1_HUMAN"}), not an accession, and
#' TTD's dedicated UniProt file (\code{P2-01-TTD_uniprot_all.txt}) was
#' checked and carries the same mnemonic, not an accession either - see
#' the file header note. Accession resolution instead goes through the
#' package's existing UniProt REST ID-mapping scaffold
#' (\code{.resolveGeneIds()}, \code{idTranslation.R}), keyed on
#' \code{GeneName} (an HGNC-style symbol, already parsed) rather than on
#' string-mangling the mnemonic (which would be wrong for TTD's many
#' non-human targets, e.g. \code{"ERG6_PNEC8"} - not \code{"_HUMAN"}).
#' Batches every distinct \code{GeneName} into a single UniProt ID-mapping
#' job (not one per row), and is wrapped so a UniProt outage degrades the
#' whole batch to \code{"resolution_failed"} rather than aborting
#' \code{\link{buildTtdDb}} entirely.
#'
#' Returns one row per \code{TargetID} with:
#' \itemize{
#'   \item \code{Uniprot_acc}: resolved accession, or \code{NA}.
#'   \item \code{uniprot_source}: \code{"TTD_missing"} (TTD's own
#'     \code{Uniprot} field was already \code{"NOUNIPROTAC"}/\code{NA} -
#'     TTD's own gap, not a mapping failure - confirmed live to be
#'     overwhelmingly family/pathway/element-level "targets" with no
#'     single accession by construction, not a mnemonic-resolution bug),
#'     \code{"unresolved"} (TTD had a real mnemonic/GeneName but the
#'     UniProt REST lookup found no accession for it - a mapping gap),
#'     \code{"resolved"}, or \code{"resolution_failed"} (the UniProt REST
#'     call itself errored/timed out - distinct from \code{"unresolved"}
#'     since it says nothing about actual UniProt coverage).
#' }
#' @keywords internal
.ttdResolveUniprotAcc <- function(targets, taxId = 9606L) {
    ttdMissing <- is.na(targets$Uniprot) | targets$Uniprot == "NOUNIPROTAC"
    acc <- rep(NA_character_, nrow(targets))
    src <- ifelse(ttdMissing, "TTD_missing", NA_character_)

    attempt <- !ttdMissing & !is.na(targets$GeneName) & nzchar(targets$GeneName)
    candidates <- unique(targets$GeneName[attempt])
    if (length(candidates)) {
        resolved <- tryCatch(
            .resolveGeneIds(candidates, idType = "symbol", to = "uniprot", taxId = taxId),
            error = function(e) NULL
        )
        if (is.null(resolved)) {
            src[attempt] <- "resolution_failed"
        } else {
            acc[attempt] <- unname(resolved[targets$GeneName[attempt]])
            src[attempt] <- ifelse(is.na(acc[attempt]), "unresolved", "resolved")
        }
    }
    ## Had a real TTD mnemonic, but no usable GeneName to look up at all.
    src[!ttdMissing & !attempt] <- "unresolved"
    data.frame(TargetID = targets$TargetID, Uniprot_acc = acc,
              uniprot_source = src, stringsAsFactors = FALSE)
}

#' small-molecule/antibody/antisense/other classifier
#'
#' Primary signal is TTD's own \code{DrugType} (\code{DRUGTYPE} field,
#' see \code{\link{.ttdParseDrugType}}), normalized from its ~80 raw
#' strings down to the 4 requested buckets by keyword match, in this
#' priority order (checked in this sequence so e.g. \code{"Antibody-drug
#' conjugate"} lands in \code{"antibody"} rather than being missed):
#' \enumerate{
#'   \item \code{DrugType} matching \code{"antisense"}/\code{"sirna"}/
#'     \code{"rnai"}/\code{"short hairpin"}/\code{"mrna"}/
#'     \code{"aptamer"}/\code{"oligonucleotide"}/\code{"nucleotide"}/
#'     \code{"nucleic acid"}/\code{"mirna"} -> \code{"antisense"};
#'   \item \code{DrugType} matching \code{"antibody"}/
#'     \code{"immunoglobulin"}/\code{"nanobody"}/\code{"probody"} ->
#'     \code{"antibody"};
#'   \item \code{DrugType} starting with \code{"small molec"} (covers
#'     \code{"Small molecular drug"}, \code{"Small molecule"}, \code{"Small
#'     molecule immunotherapy"}) -> \code{"small molecule"};
#'   \item any other non-missing \code{DrugType} (e.g. \code{"Vaccine"},
#'     \code{"CAR T Cell Therapy"}, \code{"Peptide"}, \code{"Carbohydrates"},
#'     \code{"Recombinant protein"}) -> \code{"other"}.
#' }
#' For the ~30\% of drugs with no \code{DrugType} at all, falls back to a
#' heuristic: a non-missing, non-blank \code{Smiles} value ->
#' \code{"small molecule"} (\code{Smiles} is only ever populated from
#' TTD's P1-02 \code{DRUGSMIL} field for a genuine discrete small-molecule
#' structure); else drug name/MOA text matching an antibody-nomenclature
#' suffix (\code{"...mab"}) or \code{"antibody"}/\code{"immunoglobulin"}
#' -> \code{"antibody"}; else matching \code{"antisense"}/\code{"sirna"}/
#' \code{"antagomir"}/\code{"oligonucleotide"} -> \code{"antisense"};
#' else -> \code{"other"}. Neither path is authoritative for edge cases -
#' e.g. a handful of TTD records assign a SMILES string to what is
#' actually a peptide/fusion-protein drug (confirmed live for one
#' example), which the \code{DrugType}-driven classification defers to
#' rather than second-guessing.
#' @keywords internal
.ttdMoleculeType <- function(drugType, smiles, drugName, moa) {
    dt <- trimws(drugType)
    hasType <- !is.na(dt) & nzchar(dt)

    bucket <- rep(NA_character_, length(dt))
    isAntisenseType <- grepl(paste0("antisense|sirna|rnai\\b|short hairpin|",
                                    "\\bmrna\\b|aptamer|oligonucleotide|",
                                    "nucleotide|nucleic acid|\\bmirna\\b"),
                             dt, ignore.case = TRUE)
    isAntibodyType <- grepl("antibody|immunoglobulin|nanobody|probody",
                            dt, ignore.case = TRUE)
    isSmallMolType <- grepl("^small molec", dt, ignore.case = TRUE)
    bucket[hasType] <- "other"
    bucket[hasType & isSmallMolType]  <- "small molecule"
    bucket[hasType & isAntibodyType]  <- "antibody"
    bucket[hasType & isAntisenseType] <- "antisense"

    fallback <- !hasType
    if (any(fallback)) {
        hasSmiles <- !is.na(smiles[fallback]) & nzchar(trimws(smiles[fallback]))
        text <- paste(drugName[fallback], moa[fallback])
        isAb <- grepl("mab\\b|antibody|immunoglobulin", text, ignore.case = TRUE)
        isAs <- grepl("antisense|sirna|antagomir|oligonucleotide", text, ignore.case = TRUE)
        bucket[fallback] <- ifelse(hasSmiles, "small molecule",
                            ifelse(isAb, "antibody",
                            ifelse(isAs, "antisense", "other")))
    }
    bucket
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
#' @return character(1) local file path to the SQLite database. The
#'   database additionally carries a \code{ttd_release} table (columns
#'   \code{file}, \code{version}, \code{date}) recording every downloaded
#'   TTD file's own embedded release stamp separately, since TTD does not
#'   reliably move these together across files - see the file header
#'   note. \code{\link{ttdTargetAnnot}} attaches this as
#'   \code{attr(result, "ttd_release")}.
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
    drugType    <- .ttdParseDrugType(paths$drugs)
    drugXref    <- .ttdParseDrugXref(paths$crossmatch)
    indications <- .ttdPackIndications(.ttdParseDrugIndications(paths$drugDisease))
    uniprotAcc  <- .ttdResolveUniprotAcc(targets)
    mapping     <- as.data.frame(readxl::read_excel(paths$mapping))

    interactions <- merge(mapping, targets, by = "TargetID", all.x = TRUE)
    interactions <- merge(interactions, uniprotAcc, by = "TargetID", all.x = TRUE)
    interactions <- merge(interactions, drugNames, by = "DrugID", all.x = TRUE)
    interactions <- merge(interactions, drugSmiles, by = "DrugID", all.x = TRUE)
    interactions <- merge(interactions, drugType, by = "DrugID", all.x = TRUE)
    interactions <- merge(interactions, drugXref, by = "DrugID", all.x = TRUE)
    interactions <- merge(interactions, indications, by = "DrugID", all.x = TRUE)

    ## Data-quality normalization: TTD's own literal "." MOA sentinel ->
    ## NA (confirmed live: 2253/45453 mapping rows; MOA is a plain
    ## character column here, not a factor, so no factor-level remap
    ## is needed).
    interactions$MOA[interactions$MOA == "."] <- NA_character_

    interactions$molecule_type <- .ttdMoleculeType(
        interactions$DrugType, interactions$Smiles, interactions$DrugName, interactions$MOA)

    ## Existing 9 columns first, unchanged, then the new ones - additive,
    ## same one-row-per-drug-target-edge grain as before (every new
    ## lookup below is unique-keyed on TargetID or DrugID, so none of
    ## these merges can multiply rows).
    interactions <- interactions[, c(
        "TargetID", "GeneName", "Uniprot", "TargetType",
        "DrugID", "DrugName", "Smiles", "Highest_status", "MOA",
        "Uniprot_acc", "uniprot_source", "molecule_type",
        "PubChem_CID", "PubChem_SID", "CAS", "ChEBI_ID", "Indication")]

    release <- do.call(rbind, lapply(names(paths), function(nm) {
        r <- .ttdRelease(paths[[nm]])
        data.frame(file = nm, version = r$version, date = r$date,
                  stringsAsFactors = FALSE)
    }))

    tmpDb <- tempfile(fileext = ".db")
    con <- dbConnect(SQLite(), tmpDb)
    dbWriteTable(con, "ttd_interactions", interactions, overwrite = TRUE)
    dbWriteTable(con, "ttd_release", release, overwrite = TRUE)
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
#' exact-match. Note: TTD's \code{"uniprot"} field (and the \code{col}
#' used for \code{idType = "uniprot"} lookups here) is the UniProt
#' *mnemonic entry name* (e.g. \code{"FGFR1_HUMAN"}), NOT the accession
#' number (e.g. \code{"P11362"}) used by the package's ChEMBL-SQL side
#' (\code{\link{drugTargetAnnot}}) - a real cross-source naming
#' inconsistency confirmed to persist even in TTD's own dedicated UniProt
#' file, not a bug; the separate \code{Uniprot_acc} column carries a
#' resolved accession where available (see \code{\link{buildTtdDb}} for
#' how, and its \code{uniprot_source} sibling column for why a given row
#' may still be \code{NA}), so query by \code{idType = "uniprot"} still
#' needs the mnemonic, but downstream joins to other sources can use
#' \code{Uniprot_acc} directly.
#'
#' @param queryBy list with components \code{molType}, \code{idType},
#'   \code{ids} (character vector).
#' @param ttdDbPath character(1) path to the TTD SQLite, e.g. from
#'   \code{\link{buildTtdDb}}.
#' @param fields \code{"core"} (default) or \code{"all"} - both return
#'   every \code{ttd_interactions} column, since (unlike the REST-backed
#'   sources) there is no larger raw payload to opt into - or a character
#'   vector of column names to keep (\code{QueryIDs} is always retained).
#'   See \code{\link{listDrugTargetFields}}.
#' @return A \code{data.frame} with columns \code{QueryIDs},
#'   \code{TargetID}, \code{GeneName}, \code{Uniprot}, \code{TargetType},
#'   \code{DrugID}, \code{DrugName}, \code{Smiles}, \code{Highest_status},
#'   \code{MOA}, \code{Uniprot_acc}, \code{uniprot_source},
#'   \code{molecule_type}, \code{PubChem_CID}, \code{PubChem_SID},
#'   \code{CAS}, \code{ChEBI_ID}, \code{Indication} (with the default
#'   \code{fields = "core"}), or a subset when \code{fields} requests
#'   specific columns. Also carries \code{attr(result, "ttd_release")}, a
#'   data.frame of every source TTD file's own \code{(file, version,
#'   date)} release stamp (see \code{\link{buildTtdDb}}); \code{NULL} for
#'   a \code{ttdDbPath} built before this attribute was added.
#' @examples
#' \donttest{
#'   dbPath <- buildTtdDb()
#'   ttdTargetAnnot(list(molType = "protein", idType = "symbol",
#'                        ids = c("FGFR1", "IL1B")), dbPath)
#'   ttdTargetAnnot(list(molType = "cmp", idType = "name",
#'                        ids = "Pemigatinib"), dbPath)
#' }
#' @seealso \code{\link{buildTtdDb}}, \code{\link{drugTargetAnnot}},
#'   \code{\link{listDrugTargetFields}}
#' @export
ttdTargetAnnot <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                           ttdDbPath, fields = "core") {
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
    ## NULL (not an error) for a ttdDbPath built before ttd_release
    ## existed, so older cached databases still work.
    release <- tryCatch(dbGetQuery(con, "SELECT * FROM ttd_release"),
                        error = function(e) NULL)
    out <- .dtiSelectFields(out, fields, .dtiTtdAllCols)
    attr(out, "ttd_release") <- release
    out
}

#' Documented column list for \code{listDrugTargetFields("ttd")}
#'
#' Unlike the REST-backed sources' \code{fields = "all"} (ChEMBL,
#' PubChem, DGIdb, Open Targets), \code{ttd_interactions} is a single
#' flat local SQLite table built entirely by \code{\link{buildTtdDb}}
#' (see there), so this is an exact list, not a best-effort one, and
#' \code{fields = "core"} and \code{fields = "all"} are equivalent for
#' \code{\link{ttdTargetAnnot}}.
#' @keywords internal
.dtiTtdAllCols <- c("QueryIDs", "TargetID", "GeneName", "Uniprot", "TargetType",
                    "DrugID", "DrugName", "Smiles", "Highest_status", "MOA",
                    "Uniprot_acc", "uniprot_source", "molecule_type",
                    "PubChem_CID", "PubChem_SID", "CAS", "ChEBI_ID", "Indication")
