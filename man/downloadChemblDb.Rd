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
## The ChEMBL release is large (tens of GB compressed), so it is
## downloaded into the local cache once and reused from there:
rappdirs::user_cache_dir("drugTargetInteractions")

## A small sample of the same database ships with the package, and is
## what the examples of the functions that query it use, so trying those
## out needs no download at all:
basename(system.file("extdata", "chembl_sample.db",
                     package = "drugTargetInteractions"))

\donttest{
	## The full release. Pass a version to pin one, or leave it to
	## resolve the latest live.
	downloadChemblDb(27)
}
}
