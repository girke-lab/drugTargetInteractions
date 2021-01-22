\name{downloadChemblDb}
\alias{downloadChemblDb}
\title{
downloadChemblDb
}
\description{
Download ChEMBL sqlite db for use by several other functions in the package.
}
\usage{

downloadChemblDb(version,rerun=TRUE,config=genConfig())
}
\arguments{
  \item{version}{
	The ChEMBL version to download.
	}
  \item{rerun}{
		If TRUE, the file will be downloaded, otherwise do nothing.
	}
  \item{config}{
		The configuration object. This gives the location to put the downloaded chembl db.
	}
}
\value{
 No return value.
}
\author{
Kevin Horan
}


\examples{
%% this will download a 30GB file
	\donttest{
		downloadChemblDb(27)
	}
}
