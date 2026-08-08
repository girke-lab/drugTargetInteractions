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

.writeFixture <- function(result, name) {
    saveRDS(list(result = result, cachedDate = as.character(Sys.Date())),
            file.path(fixtureDir, name))
    message("wrote ", name)
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
.writeFixture(getEnsemblParalogs("FGFR1"), "ensembl_paralogs_fgfr1.rds")
genes8 <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
.writeFixture(getEnsemblParalogs(genes8), "ensembl_paralogs_batch8genes.rds")
