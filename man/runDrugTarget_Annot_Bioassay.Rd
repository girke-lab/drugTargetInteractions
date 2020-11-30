\name{runDrugTarget_Annot_Bioassay}
\alias{runDrugTarget_Annot_Bioassay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
runDrugTarget_Annot_Bioassay
}
\description{
Meta-function to obtain in one step both drug-target annotation and bioassay data.
}
\usage{
runDrugTarget_Annot_Bioassay(res_list, up_col_id="ID", ens_gene_id,  cmpid_file=file.path(config$resultsPath,"cmp_ids.rds") ,config=genConfig(), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res_list}{
      Object obtained from \code{getUniprotIDs} function. 
}
  \item{up_col_id}{
      Column name in \code{data.frames} of \code{res_list} containing uniprot IDs, usually one of: 'ID', 'ID_up_sp' or 'ID_up_sp_tr'.
}
  \item{ens_gene_id}{
      Named character vector with ENSEMBL gene IDs in name slot and gene symbols or other ID type in value slot
}
  \item{cmpid_file}{
     Path to CMP ID mapping file, often named \code{cmp_ids.rds}.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
  \item{\dots}{
     Slot to pass on additional arguments. 
}
}
\value{
List with two components each containing a \code{data.frame}. The first one (\code{Annotation}) 
contains drug-target annotation data, and the second one (\code{Bioassay}) contains drug-target
bioassay data. 
}
\references{
References to be added...
}
\author{
Thomas Girke
}

\seealso{
See also: \code{drugTargetAnnot} and \code{drugTargetBioactivity}
}
\examples{

\dontrun{

## (1) Translate gene symbols to ENSEMBL gene IDs
ensembl_gene_id <- c("ENSG00000001626", "ENSG00000168748")
idMap <- getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=ensembl_gene_id, idtype="ENSEMBL_GENE_ID")
ens_gene_id <- idMap$ens_gene_id


## (2a) Retrieve UniProt IDs with both IDMs and SSNN paralogs
queryBy <- list(molType="gene", idType="ensembl_gene_id", ids=names(ens_gene_id))
res_list <- getParalogs(queryBy)

## (2b) Retrieve UniProt IDs with both IDMs and SSNN paralogs
res_list90 <- getUniprotIDs(taxId=9606, kt="ENSEMBL", keys=names(ens_gene_id), seq_cluster="UNIREF90", chunksize=1)

## (3) Obtain Drug-Target Annotation and Bioassay Data 
drug_target_list <- runDrugTarget_Annot_Bioassay(res_list=res_list, up_col_id="ID_up_sp", ens_gene_id ) 

}
}
\keyword{utilities}
