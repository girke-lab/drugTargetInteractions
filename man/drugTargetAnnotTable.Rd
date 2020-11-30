\name{drugTargetAnnotTable}
\alias{drugTargetAnnotTable}
\title{
drugTargetAnnotTable
}
\description{
	Generates a drug target annotation TSV file. This file
	includes target information from ChEMBL, drugbank, pubchem, and chembi.

	This function requires the ID mapping files "src1src2.txt.gz", "src1src22.txt.gz", and
	"src1src7.txt.gz" to exist in a directory called "downloads" before being run. These
	can be generated with the \code{\link{downloadUniChem}} function.
}
\usage{
drugTargetAnnotTable(outfile, rerun=TRUE,config=genConfig()) 
}
\arguments{
  \item{outfile}{
	  The name of the output file to write the results to.
}
  \item{rerun}{
	  If true, download and generate output file. Otherwise do nothing.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Writes output file to \code{outfile}.
}
\author{
	Thomas Girke
}


\seealso{
	\code{\link{downloadUniChem}}
}
\examples{
	\dontrun{
		drugTargetAnnotTable(outfile="drugTargetAnnot.xls", rerun=FALSE)
	}	

}
