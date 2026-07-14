\name{downloadChemblDb}
\alias{downloadChemblDb}
\title{
downloadChemblDb
}
\description{
Download ChEMBL sqlite db for use by several other functions in the package.
}
\usage{

downloadChemblDb(version=.chemblLatestVersion(),rerun=TRUE,config=genConfig())
}
\arguments{
  \item{version}{
	The ChEMBL version to download. Defaults to the latest release,
	resolved live from ChEMBL's REST \code{status} endpoint (see
	\code{\link{getChemblMolecule}} and related REST accessors); pass an
	explicit integer to pin to an older release instead.
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
