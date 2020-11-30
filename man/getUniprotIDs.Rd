\name{getUniprotIDs}
\alias{getUniprotIDs}
\title{
Retrieve UniProt IDs via ID and Cluster Mappings
}
\description{
The following returns for a set of query IDs (e.g. Ensembl gene IDs) the
corresponding UniProt IDs based on two independent approaches: ID mappings
(IDMs) and sequence similarity nearest neighbors (SSNNs) using UNIREF clusters.
Note, the 'keys' or query IDs (e.g. ENSEMBL genes) can only be reliably maintained
in the SSNN results when 'chunksize=1' since batch queries for protein clusters
with 'UnitProt.ws' will often drop the query IDs. To 
address this, the query result contains an extra 'QueryID' column when 'chunksize=1',
but not when it is set to a different value than 1.

The \code{\link{getParalogs}} function is similar but it uses biomaRt's paralogs instead of UNIREF clusters.
}
\usage{
getUniprotIDs(taxId = 9606, kt = "ENSEMBL", keys, seq_cluster = "UNIREF90", chunksize=20)
}
\arguments{
  \item{taxId}{
	  An NCBI taxonomy ID
}
  \item{kt}{
	  Should be either  "ENSEMBL" or "UNIPROTKB".
}
  \item{keys}{
	  Query IDs.
}
  \item{seq_cluster}{
%%     ~~Describe \code{seq_cluster} here~~
	Which cluster to use. Should be one of  'UNIREF100', 'UNIREF90', 'UNIREF50'.
}
  \item{chunksize}{
	  Queries are done in batches, this parameter sets the size of each batch.
}
}
\value{
	Returns a list of data. 
}
\author{
	Thomas Girke
}


\seealso{
 \code{\link{getParalogs}} 
 \code{\link{UniProt.ws}}
}
\examples{
%%  requires network 
	\dontrun{ 
		 keys <- c("ENSG00000145700", "ENSG00000135441", "ENSG00000120071", "ENSG00000120088", "ENSG00000185829", "ENSG00000185829", "ENSG00000185829", "ENSG00000238083", "ENSG00000012061", "ENSG00000104856", "ENSG00000104936", "ENSG00000117877", "ENSG00000130202", "ENSG00000130202", "ENSG00000142252", "ENSG00000189114", "ENSG00000234906") 
		 res_list100 <- getUniprotIDs(taxId=9606, kt="ENSEMBL", keys=keys, seq_cluster="UNIREF100") 
	}

}
