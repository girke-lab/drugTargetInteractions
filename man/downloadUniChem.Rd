\name{downloadUniChem}
\alias{downloadUniChem}
\title{
downloadUniChem
}
\description{
	Downloads UniChem compound ID mappings from https://www.ebi.ac.uk/unichem/ucquery/listSources.
	Mappings are downloaded for DrugBank, PubChem, and ChEBI.
}
\usage{
downloadUniChem(rerun=TRUE, config=genConfig())
}
\arguments{
  \item{rerun}{
	  If true, downloads the files, else does nothing.
	}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Generates the following output files: "src1src2.txt.gz", "src1src22.txt.gz", and
	"src1src7.txt.gz". These correspond to mappings from ChEMBL to DrugBank, PubChem, and 
	ChEBI, respectivly. 
}
\references{
https://www.ebi.ac.uk/unichem/ucquery/listSources
}
\author{
	Thomas Girke
}


\seealso{

	\code{\link{drugTargetAnnotTable}}
}
\examples{
	\dontrun{
		downloadUniChem(rerun=TRUE)
	}
}
