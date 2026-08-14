## =====================================================================
##  generate_live_fixtures.R
##  Maintainer-run script that (re)generates the shipped fallback
##  fixtures under inst/extdata/fixtures/, used by .dtiLiveOrCached()
##  (R/liveFallback.R) when a live API call fails in a test or vignette
##  chunk.
##
##  Run manually, NOT as part of any automated build:
##    - whenever a source's schema is suspected to have changed
##    - periodically, to keep the shipped example/fixture data reasonably
##      current
##
##  Each fixture is an .rds holding list(result = <object>, cachedDate =
##  "YYYY-MM-DD") - metadata lives alongside the object rather than as
##  attributes on it, so it survives independently of whatever the live
##  function itself returns.
##
##  Run from the package root with the dev source loaded, e.g.:
##    devtools::load_all(); source("inst/scripts/generate_live_fixtures.R")
## =====================================================================

fixtureDir <- "inst/extdata/fixtures"
if (!dir.exists(fixtureDir)) dir.create(fixtureDir, recursive = TRUE)

## Refuse to write an empty fixture. The transport helpers degrade to
## NULL/empty on an upstream failure rather than erroring, so a run
## started during an outage will otherwise happily overwrite a good
## fixture with nothing - which is exactly how
## ensembl_paralogs_fgfr1.rds came to ship with 0 rows on 2026-08-08,
## silently defeating the fallback it was meant to provide. An empty
## result here is always a failed fetch, never a real answer: every
## fixture below is a query chosen precisely because it returns data.
## `check` guards against *partial* capture, which an emptiness test
## cannot catch and which is the more insidious form: a batch query
## whose per-item requests fail individually still returns a non-empty
## frame, just one quietly missing items. That is how
## ensembl_paralogs_batch8genes.rds came to ship covering 6 of its 8
## genes. Pass a predicate asserting what a *complete* capture looks
## like for that specific fixture.
.writeFixture <- function(result, name, check = NULL) {
    n <- if (is.data.frame(result)) nrow(result)
         else if (is.list(result)) sum(vapply(result, NROW, integer(1)))
         else NROW(result)
    if (is.null(result) || n == 0L)
        stop("refusing to write empty fixture '", name, "' - the live call ",
             "returned nothing, which means the source was failing. Re-run ",
             "when it is healthy; the existing fixture is left untouched.")
    if (!is.null(check) && !isTRUE(check(result)))
        stop("refusing to write incomplete fixture '", name, "' - the live ",
             "call returned ", n, " row(s) but failed its completeness ",
             "check, so some part of the query silently did not come back. ",
             "Re-run when the source is healthy; the existing fixture is ",
             "left untouched.")
    saveRDS(list(result = result, cachedDate = as.character(Sys.Date())),
            file.path(fixtureDir, name))
    message("wrote ", name, " (", n, " rows)")
}

## --- ChEMBL: vignette 'fields_all' chunk (CHEMBL25 = aspirin) --------------
chembl_all <- getChemblDrugTarget(
    list(molType = "cmp", idType = "chembl_id", ids = "CHEMBL25"),
    fields = "all")
.writeFixture(chembl_all, "chembl_fields_all_chembl25.rds")

## --- Legacy biomaRt getParalogs(): vignette 'legacy_run_everything' -------
idMap <- getSymEnsUp(EnsDb = "EnsDb.Hsapiens.v86",
                     ids = c("CA7", "CFTR"), idtype = "GENE_NAME")
queryBy <- list(molType = "gene", idType = "ensembl_gene_id",
               ids = names(idMap$ens_gene_id))
paralogsLegacy <- getParalogs(queryBy)
.writeFixture(paralogsLegacy, "ensembl_legacy_paralogs_ca7_cftr.rds")

## --- Ensembl REST: test-idTranslation.R's 5 tests -------------------------
.writeFixture(getEnsemblParalogs("NLRP3"), "ensembl_paralogs_nlrp3.rds")
.writeFixture(getEnsemblOrthologs("NLRP3", targetSpecies = "mouse"),
             "ensembl_orthologs_nlrp3_mouse.rds")
.writeFixture(getEnsemblParalogs("NLRP3", condensed = TRUE),
             "ensembl_paralogs_nlrp3_condensed.rds")
.writeFixture(getEnsemblParalogs("FGFR1"), "ensembl_paralogs_fgfr1.rds",
             check = function(r) nrow(r) >= 50L)
genes8 <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
## Every one of the 8 must be represented, or the batch was only
## partially fetched - see .writeFixture()'s note.
.writeFixture(getEnsemblParalogs(genes8), "ensembl_paralogs_batch8genes.rds",
             check = function(r) setequal(unique(r$query_gene), genes8))
