\name{cmpIdMapping}
\alias{cmpIdMapping}
\title{
cmpIdMapping
}
\description{
	Function to generate compound ID mappings UniChem. 

	This function requires the ID mapping files "src1src2.txt.gz", "src1src22.txt.gz", and
	"src1src7.txt.gz" to exist in a directory called "downloads" before being run. These
	can be generated with the \code{\link{downloadUniChem}} function.

	It will do some processing on these files and output an RDS file
	at \code{outfile}. This file can then be used in other functions,
	such as \code{\link{drugTargetAnnot}}.

}
\usage{
cmpIdMapping(outfile=file.path(config$resultsPath,"cmp_ids.rds"), rerun=TRUE,config=genConfig()) 
}
\arguments{
  \item{outfile}{
	  Path to output file.
}
  \item{rerun}{
	  If true, runs processing, otherwise does nothing.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Generates an RDS file at \code{outfile}.
}
\author{
	Thomas Girke
}


\seealso{
	\code{\link{downloadUniChem}}
}
\examples{
	\dontrun{
		 cmpIdMapping(outfile="cmp_ids.rds", rerun=TRUE)
		 cmp_ids <- readRDS("cmp_ids.rds")
	}
}
