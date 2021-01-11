\name{genConfig}
\alias{genConfig}
\title{
genConfig
}
\description{
Create a default configuration object. 
}
\usage{

genConfig( chemblDbPath = "chembldb.db", downloadPath = "downloads", resultsPath = "results" )
}
\arguments{
  \item{chemblDbPath}{
		Path or filename of ChEMBL SQLite db file.
	}
  \item{downloadPath}{
		The name of a directory to put downloaded files in.	
	}
  \item{resultsPath}{
		The name of a directory to put output files in.	
	}
}
\value{
	A config object that can be passed to ther functions.
}
\author{
Kevin Horan
}


\examples{
	config = genConfig()
}
