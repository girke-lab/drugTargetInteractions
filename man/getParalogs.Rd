\name{getParalogs}
\alias{getParalogs}
\title{
	getParalogs
}
\description{
Using \code{biomaRt}, obtain for query genes the corresponding UniProt IDs as well
as paralogs.  Query genes can be Gene Names or ENSEMBL Gene IDs from 
H sapiens. The result is similar to IDMs and SSNNs from \code{\link{getUniprotIDs}}
function, but instead of UNIREF clusters, biomaRt's paralogs are used to
obtain SSNNs. 

}
\usage{
getParalogs(queryBy)
}
\arguments{
  \item{queryBy}{
	  A list defining the query, as described in  \code{\link{queryBy}}.
}
}
\value{
	Returns a list with the paralogs for the given genes.
}
\author{
	Thomas Girke
}


\seealso{
	\code{\link{getUniprotIDs}}
	\code{\link{queryBy}}
}
\examples{
	\dontrun{
		 queryBy <- list(molType="gene", idType="external_gene_name", ids=c("ZPBP", "MAPK1", "EGFR"))
		 result <- getParalogs(queryBy)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
