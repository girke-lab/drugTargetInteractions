\name{getSymEnsUp}
\alias{getSymEnsUp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene to Protein ID Mappings
}
\description{
The \code{getSymEnsUp} function returns for a query of gene or protein IDs a mapping
table containing: ENSEMBL Gene IDs, Gene Names/Symbols, UniProt IDs and ENSEMBL 
Protein IDs. Subsequent slots contain the corresponding named character vectors.
Internally, the function uses the \code{ensembldb} package. 
}
\usage{
getSymEnsUp(EnsDb = "EnsDb.Hsapiens.v86", ids, idtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EnsDb}{
\code{EnsDb} instance of ensembldb package
}
  \item{ids}{
Character vector with IDs matching the type specified under \code{idtype}
}
  \item{idtype}{
Character vector of length one containing one of: \code{GENE_NAME}, \code{ENSEMBL_GENE_ID} or \code{UNIPROT_ID}
}
}
\value{
List object with following components:
    \item{idDF}{ID mapping \code{data.frame}}
    \item{ens_gene_id}{named character vector}
    \item{up_ens_id}{named character vector}
    \item{up_gene_id}{named character vector}
}
\author{
Thomas Girke
}

\examples{
	 gene_name <- c("CA7", "CFTR")
	 getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=gene_name, idtype="GENE_NAME")
}
\keyword{ utilities }
