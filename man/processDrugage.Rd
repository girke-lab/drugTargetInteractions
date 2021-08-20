\name{processDrugage}
\alias{processDrugage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
processDrugage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Download Drug Age data from genomics.senescence.info/drugs. Process the data
and write it out as a TSV spreadsheet.
}
\usage{
processDrugage(drugagefile=file.path(config$resultsPath,"drugage_id_mapping.xls"), redownloaddrugage=TRUE,config=genConfig()) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drugagefile}{
	  The name of the output file.
}
  \item{redownloaddrugage}{
	  If true, download the data file. Otherwise assume the file is already downloaded.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Output is written to \code{drugagefile}.
}
\author{
	Thomas Girke
}


\examples{

	tryCatch({
		config = genConfig(chemblDbPath= 
									system.file("extdata", "chembl_sample.db", package="drugTargetInteractions"))
		processDrugage("druage_id_mapping.xls",TRUE,config)
	},
	error=function(e){
			message("Failed to run processDrugage(), please try again later")
		}
	)
}
