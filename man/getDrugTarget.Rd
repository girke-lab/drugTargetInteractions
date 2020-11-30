\name{getDrugTarget}
\alias{getDrugTarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getDrugTarget
}
\description{
	This function allows you to query a subset of the data 
	fetched by \code{\link{drugTargetAnnotTable}}.
}
\usage{
getDrugTarget(dt_file=file.path(config$resultsPath,"drugTargetAnnot.xls"), queryBy=list(molType=NULL, idType=NULL, ids=NULL), 
								  id_mapping=c(chembl="chembl_id", pubchem="PubChem_ID", uniprot="UniProt_ID"), columns,config=genConfig()) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dt_file}{
	  The drug target annotation file. This can be generated with \code{\link{drugTargetAnnotTable}}.
}
  \item{queryBy}{
	  A list defining the query, as described in  \code{\link{queryBy}}.
}
  \item{id_mapping}{
	  A list providing the id columns for ChEMBL, PubChem, and UniProt. It should contain the 
	  fields "chembl", "pubchem", and "uniprot", each wit the column name of the respective id number
	  in the drug target annotation file.  See default value above for an example.

}
  \item{columns}{
	  A list of column indexes to select as a subset of the final result set.
}
	\item{config}{ General configuration. See \code{\link{genConfig}}.  }
}
\value{
	Returns the query result as a data frame. 
}
\author{
	Thomas Girke
}

\seealso{
	\code{\link{drugTargetAnnotTable}}
}
\examples{
	\dontrun{
		 id_mapping <- c(chembl="chembl_id", pubchem="PubChem_ID", uniprot="UniProt_ID", drugbank="DrugBank_ID")

		 queryBy <- list(molType="cmp", idType="chembl", ids=c("CHEMBL25", "CHEMBL1742471"))
		 getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52)) 

		 queryBy <- list(molType="cmp", idType="pubchem", ids=c("2244", "65869", "2244"))
		 getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52)) 

		 queryBy <- list(molType="protein", idType="uniprot", ids=c("P43166", "P00915", "P43166"))
		 getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52))


	}
}
