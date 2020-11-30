\name{transformTTD}
\alias{transformTTD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transformTTD
}
\description{
Integration with Therapeutic Target Database (TTD).
This function downloads a data file from idrblab.org and
returns it as a data frame.
}
\usage{
transformTTD(ttdfile=file.path(config$downloadPath,"TTD_IDs.txt"), redownloadTTD=TRUE,config=genConfig()) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttdfile}{
	  The name of the output file to write the downloaded file to.
}
  \item{redownloadTTD}{
	  If true, data file will be downloaded again. If false, we assume the file
	  already exists at \code{ttdfile}.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Returns a data frame with TTD data in it.
}
\author{
	Thomas Girke
}


\examples{
	\dontrun{
		ttd = transformTTD("ttd.txt",TRUE)
	}

}
