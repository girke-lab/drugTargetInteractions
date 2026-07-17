## Network-guarded tests for the direct-REST ID-translation layer
## (UniProt ID mapping + Ensembl homology). Defined locally (not shared
## from test-apiAccess.R) so this file runs standalone too.

skip_if_offline_dti <- function() {
    testthat::skip_on_cran()
    if (!.dtiHasInternet())
        testthat::skip("No internet / API unreachable")
}

test_that("getUniprotMapping resolves gene symbols to Swiss-Prot accessions (human), excludes unmapped IDs", {
    skip_if_offline_dti()
    m <- getUniprotMapping(c("FGFR1", "KLB", "NOTAREALGENEXYZ"),
                           from = "Gene_Name", to = "UniProtKB-Swiss-Prot")
    expect_s3_class(m, "data.frame")
    expect_identical(names(m), c("From", "To"))
    expect_identical(m$To[m$From == "FGFR1"], "P11362")
    expect_identical(m$To[m$From == "KLB"], "Q86Z14")
    expect_false("NOTAREALGENEXYZ" %in% m$From)  # no per-ID failure row, by design
})

test_that("getUniprotMapping's taxId filter actually restricts to one organism", {
    skip_if_offline_dti()
    ## OLFM1 has Swiss-Prot entries in multiple species; taxId=9606 must
    ## return exactly the human one (this is the exact scenario that was
    ## silently broken in UniProt.ws's mapUniProt() until v2.49.1).
    m <- getUniprotMapping("OLFM1", from = "Gene_Name",
                           to = "UniProtKB-Swiss-Prot", taxId = 9606L)
    expect_equal(nrow(m), 1L)
    expect_identical(m$To, "Q99784")
})

test_that("getUniprotMapping reverse direction (accession -> Ensembl) carries a version suffix", {
    skip_if_offline_dti()
    m <- getUniprotMapping("P11362", from = "UniProtKB_AC-ID", to = "Ensembl")
    expect_true(grepl("^ENSG00000077782\\.", m$To))
})

test_that("getUniprotMapping paginates through more than one results page", {
    skip_if_offline_dti()
    genes <- c("FGFR1", "FGFR2", "FGFR3", "FGFR4", "EGFR", "ERBB2", "MET",
              "KIT", "PDGFRA", "PDGFRB", "ALK", "ROS1", "RET", "NTRK1",
              "NTRK2", "NTRK3")
    m <- getUniprotMapping(genes, from = "Gene_Name", to = "UniProtKB", taxId = NULL)
    expect_gt(nrow(m), 500L)  # default page size; confirms pagination followed Link headers
    expect_setequal(unique(m$From), genes)
})

test_that("getUniprotMapping rejects empty ids", {
    expect_error(getUniprotMapping(character(0), from = "Gene_Name", to = "Ensembl"))
})

test_that(".dtiUniprotParseResults handles both plain-string and nested-object 'to' shapes, no network", {
    results <- list(
        list(from = "FGFR1", to = "P11362"),
        list(from = "KLB", to = list(primaryAccession = "Q86Z14"))
    )
    df <- .dtiUniprotParseResults(results)
    expect_equal(nrow(df), 2L)
    expect_identical(df$To, c("P11362", "Q86Z14"))
    expect_identical(.dtiUniprotParseResults(list()), .dtiUniprotParseResults(NULL))
})

test_that("getEnsemblParalogs returns tidy within-human paralogs for NLRP3", {
    skip_if_offline_dti()
    p <- getEnsemblParalogs("NLRP3")
    expect_s3_class(p, "data.frame")
    expect_identical(names(p), c("query_gene", "homolog_id", "homolog_protein_id",
                                 "homolog_species", "type", "taxonomy_level",
                                 "perc_id", "perc_pos"))
    expect_true(nrow(p) >= 1L)
    expect_true(all(p$query_gene == "NLRP3"))
    expect_true(all(p$homolog_species == "homo_sapiens"))
    expect_true(all(!is.na(p$perc_id)))  # cigar_line=0 keeps perc_id, unlike condensed
})

test_that("getEnsemblOrthologs restricted to mouse returns NLRP3's one2one ortholog", {
    skip_if_offline_dti()
    o <- getEnsemblOrthologs("NLRP3", targetSpecies = "mouse")
    expect_equal(nrow(o), 1L)
    expect_identical(o$homolog_species, "mus_musculus")
    expect_identical(o$homolog_id, "ENSMUSG00000032691")
    expect_identical(o$type, "ortholog_one2one")
})

test_that("getEnsemblParalogs(condensed=TRUE) is faster-path: no perc_id, but type/id present", {
    skip_if_offline_dti()
    p <- getEnsemblParalogs("NLRP3", condensed = TRUE)
    expect_true(nrow(p) >= 1L)
    expect_true(all(is.na(p$perc_id)))
    expect_true(all(!is.na(p$homolog_id)))
    expect_true(all(!is.na(p$type)))
})

test_that("getEnsemblParalogs handles a large gene family (FGFR1) without error", {
    skip_if_offline_dti()
    ## Regression guard for the pathological-payload finding: default
    ## query shape (cigar_line=0) must stay well clear of .dtiApiGET's
    ## timeout even for a gene with dozens of paralogs.
    p <- getEnsemblParalogs("FGFR1")
    expect_equal(nrow(p), 53L)
})

test_that("getEnsemblParalogs/getEnsemblOrthologs batch across multiple genes, skipping those with none", {
    skip_if_offline_dti()
    genes <- c("FGF21", "KLB", "FGFR1", "NLRP3", "IL1B", "TFEB", "ADIPOR1", "ADIPOR2")
    p <- getEnsemblParalogs(genes)
    expect_true(setequal(unique(p$query_gene), genes))  # all 8 have >=1 paralog live-confirmed
})
