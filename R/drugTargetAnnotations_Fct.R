#########################################
## Helper Functions for ChEMBL Queries ##
#########################################
## Author: Thomas Girke
## Last update: 05-Apr-20

genConfig<- function(
	chemblDbPath = "chembldb.db",
	downloadPath = "downloads",
	resultsPath = "results" ){

	if(!dir.exists(downloadPath))
		dir.create(downloadPath)
	if(!dir.exists(resultsPath))
		dir.create(resultsPath)
	return(list(chemblDbPath=chemblDbPath,downloadPath=downloadPath,resultsPath=resultsPath))

}

#############
## UniChem ##
#############
## UniChem CMP ID mappings from here: 
##      https://www.ebi.ac.uk/unichem/ucquery/listSources
## Note: above html table gives numbering to select proper src_id for ftp downloads, e.g. DrugBank is src2
##      ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/
## Examples:
downloadUniChem <- function(rerun=TRUE,config = genConfig()) {
    if(rerun==TRUE) {
        ## ChEMBL to DrugBank mapping in src1src2.txt.gz: 
        download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz", file.path(config$downloadPath,"src1src2.txt.gz"))
        ## ChEMBL to PubChem CID mapping in src1src22.txt.gz
        download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz", file.path(config$downloadPath,"src1src22.txt.gz"))
        ## ChEMBL to ChEBI mapping in src1src7.txt.gz
        download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src7.txt.gz", file.path(config$downloadPath,"src1src7.txt.gz"))
    }
}
downloadChemblDb <- function(version,rerun=TRUE,config=genConfig()){
	if(rerun==TRUE){
		url = paste("ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_",version,"_sqlite.tar.gz",sep="")
		tarFile = file.path(config$downloadPath,"chembl_sqlite.tar.gz")
		tempDir = tempdir()

		dbPath = paste("chembl_",version,"/chembl_",version,"_sqlite/chembl_",version,".db",sep="")
		download.file(url,tarFile)
		untar(tarFile,c(dbPath),exdir=tempDir)
		sourceDbFile = file.path(tempDir,
										 paste("chembl_",version,sep=""),
										 paste("chembl_",version,"_sqlite",sep=""),
										 paste("chembl_",version,".db",sep=""))

		file.copy(from=sourceDbFile, to = config$chemblDbPath)
		file.remove(sourceDbFile)
	}
}
## Usage:
#downloadUniChem(rerun=FALSE)

###############################################################
## Retrieve with UniProt.ws UniProt IDs using IDMs and SSNNs ## 
###############################################################
## The following returns for a set of query IDs (here Ensembl gene or UniProt IDs) 
## the corresponding UniProt IDs based on two independent approaches: ID mappings
## (IDMs) and sequence similarity nearest neighbors (SSNNs) using UNIREF
## clusters. Note: the query IDs (e.g. ENSEMBL genes) can only be reliably maintained
## in the SSNN results when 'chunksize=1' since batch queries for protein clusters
## with UnitProt.ws will often drop the query IDs such as ENSEMBL gene IDs. To 
## address this, the query result contains an extra 'QueryID' column when 'chunksize=1',
## but not when it is set to a different value than 1.
## The getParalogs function is similar but it uses biomaRt's paralogs instead of UNIREF clusters.

## Basic usage of UniProt.ws package
# library(UniProt.ws)
# up <- UniProt.ws(taxId=9606) # Attach organism here human
# ?UniProt.ws # Help on base class of package
# columns(up) # Lists available data 
# keytypes(up) # Lists available keys
# species(up) # Returns attached species
# taxId(up) # Returns attached species id
# availableUniprotSpecies()[1:4,] # Lists available species
# availableUniprotSpecies("musculus") # Pattern matching on species
# lookupUniprotSpeciesFromTaxId("3702") # Returns species name for species id
# # taxId(up) <- 10090 # Attaches mouse as species to UniProt.ws (here 'up') object

getUniprotIDs <- function(taxId=9606, kt="ENSEMBL", keys, seq_cluster="UNIREF90", chunksize=20) {
    up <- UniProt.ws(taxId) # Attach organism, here human
    
    ## Validity checks
    if(!kt %in% c("ENSEMBL", "UNIPROTKB")) stop("Argument kt needs to be one of: 'ENSEMBL', 'UNIPROTKB'")
    if(!seq_cluster %in% c("UNIREF100", "UNIREF90", "UNIREF50")) stop("Argument seq_cluster needs to be one of: 'UNIREF100', 'UNIREF90', 'UNIREF50'")

    ## Columns to return 
    df <- data.frame('ENSEMBL'=character(), 
                     'ID'=character(), 
                     'GENES'=character(), 
                     'UNIREF100'=character(), 
                     'UNIREF90'=character(), 
                     'UNIREF50'=character(), 
                     'ORGANISM'=character(), 
                     'PROTEIN-NAMES'=character(), 
                     check.names=FALSE)
    columns <- colnames(df)

    ## Arrange keys in chunk list, with list components of length <= 'chunksize'
    keys <- unique(keys)
    keylist <- split(keys, ceiling(seq_along(keys)/chunksize))  
    
    ## Empty list container
    listcontainer <- vector("list", length(keylist))
    names(listcontainer) <- names(keylist)
    res_list <- list(IDM=listcontainer, SSNN=listcontainer)

    ## Run queries in chunks
    for(i in seq_along(keylist)) {
        ## ID mappings (IDMs)
        res_ID <- try(select(up, keys=keylist[[i]], columns, kt), silent=TRUE)
        if(!inherits(res_ID, "try-error")) {
            if(chunksize==1) {
                res_ID <- data.frame(QueryID=keylist[[i]][1], res_ID, check.names=FALSE) 
                res_list[[1]][[i]] <- res_ID[,c("QueryID", columns)]
            } else {
                res_list[[1]][[i]] <- res_ID[,columns]
            }
        } else {
            if(chunksize==1) {
                res_list[[1]][[i]] <- data.frame(QueryID=character(), df, check.names=FALSE)
                res_list[[2]][[i]] <- data.frame(QueryID=character(), df, check.names=FALSE) # Required since current iteration will not enter next step (SSNN)
            } else {
                res_list[[1]][[i]] <- df
                res_list[[2]][[i]] <- df # Required since current iteration will not enter next step (SSNN)
            }
        }

        ## Sequence Similarity (SSNNs) using UNIREF*
        if(!inherits(res_ID, "try-error")) {
            unirefkeys <- as.character(na.omit(unique(res_ID[,seq_cluster])))
            unirefkeylist <- split(unirefkeys, ceiling(seq_along(unirefkeys)/chunksize))  
            unirefcontainer <- vector("list", length(unirefkeylist))
            names(unirefcontainer) <- names(unirefkeylist)
            ## SSNN queries are run in subloop since number of uniref keys is larger than in IDM query
            for(j in seq_along(unirefkeylist)) {
                tmp <- try(select(up, keys=unirefkeylist[[j]], columns, seq_cluster), silent=TRUE)
                if(!inherits(tmp, "try-error")) {
                    if(chunksize==1) {
                        unirefcontainer[[j]] <- data.frame(QueryID=keylist[[i]][1], tmp[,columns], check.names=FALSE) 
                    } else {
                        unirefcontainer[[j]] <- tmp[,columns]
                    }
                } else {
                    if(chunksize==1) {
                        unirefcontainer <- data.frame(QueryID=character(), df, check.names=FALSE)
                    } else {
                        unirefcontainer[[j]] <- df
                    }
                }
            }
            res_CL <- do.call("rbind", unirefcontainer)
            if(chunksize==1) {
                res_list[[2]][[i]] <- res_CL[,c("QueryID", columns)]
            } else {
                res_list[[2]][[i]] <- res_CL[,columns]
            }
        }
        ## Status message
        print(paste("Processed", chunksize * i, "keys of total of", length(keys), "keys."))
    }
    ## Assemble and return results
    res_list[[1]] <- do.call("rbind", res_list[[1]])
    rownames(res_list[[1]]) <- NULL
    res_list[[2]] <- do.call("rbind", res_list[[2]])
    rownames(res_list[[2]]) <- NULL
    return(res_list)
}

## Usage:
# keys <- c("ENSG00000145700", "ENSG00000135441", "ENSG00000120071", "ENSG00000120088", "ENSG00000185829", "ENSG00000185829", "ENSG00000185829", "ENSG00000238083", "ENSG00000012061", "ENSG00000104856", "ENSG00000104936", "ENSG00000117877", "ENSG00000130202", "ENSG00000130202", "ENSG00000142252", "ENSG00000189114", "ENSG00000234906") 
# res_list100 <- getUniprotIDs(taxId=9606, kt="ENSEMBL", keys=keys, seq_cluster="UNIREF100") 
# sapply(res_list100, dim, simplify=FALSE)
# sapply(names(res_list100), function(x) unique(na.omit(res_list100[[x]]$ID)))
# res_list90 <- getUniprotIDs(taxId=9606, kt="ENSEMBL", keys=keys, seq_cluster="UNIREF90") 
# sapply(res_list90, dim, simplify=FALSE)
# sapply(names(res_list90), function(x) unique(na.omit(res_list90[[x]]$ID)))

##########################################################
## Retrieve with biomaRt UniProt IDs and their Paralogs ##
##########################################################
## Using biomaRt, obtain for query genes the corresponding UniProt IDs as well
## as paralogs.  Query genes can be Gene Names or ENSEMBL Gene IDs from 
## H sapiens. The result is similar to IDMs and SSNNs from getUniprotIDs
## function, but instead of UNIREF clusters, biomaRt's paralogs are used to
## obtain SSNNs. 

getParalogs <- function(queryBy) {
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")                                                                                                            
    
    ## ID Matching (IDM) result table
    ## To list available uniprot annotation fields, run:
    # listAttributes(mart)[grep("uniprot", listAttributes(mart)[,"name"]),]
    IDMresult <- getBM(attributes = c("ensembl_gene_id",                                                                                                                        
                                    "uniprot_gn_symbol",
                                    "uniprotswissprot",
                                    "uniprotsptrembl",
                                    "description"),
                    filters = queryBy$idType,                                                                                                                          
                    values = queryBy$ids,                                                                                                                                        
                    mart = mart)                
    IDMresult <- data.frame(QueryID=as.character(IDMresult$ensembl_gene_id), 
                            ENSEMBL=as.character(IDMresult$ensembl_gene_id), 
                            GENES=as.character(IDMresult$uniprot_gn_symbol), 
                            ID_up_sp=as.character(IDMresult$uniprotswissprot), 
                            ID_up_sp_tr=as.character(IDMresult$uniprotsptrembl)) 
   
    ## Paralog (Sequence Similarity Nearest Neighbors - SSNN) result table
    ## To list available paralog annotation fields, run:
    # listAttributes(mart)[grep("paralog", listAttributes(mart)[,"name"]),]
    result <- getBM(attributes = c("external_gene_name",                                                                                                                     
                                    "ensembl_gene_id",                                                                                                                        
                                    # "hsapiens_paralog_canonical_transcript_protein",                                                                                          
                                    "hsapiens_paralog_associated_gene_name",                                                                                                  
                                    "hsapiens_paralog_ensembl_gene",                                                                                                          
                                    # "hsapiens_paralog_ensembl_peptide",
                                    "hsapiens_paralog_perc_id",
                                    "hsapiens_paralog_perc_id_r1"),                                                                                                      
                    filters = queryBy$idType,                                                                                                                          
                    values = queryBy$ids,                                                                                                                                        
                    mart = mart)                                                                                                                                            
    ## Reshape result table
    query <- result[!duplicated(result$ensembl_gene_id),]
    query <- data.frame(QueryID=as.character(query$ensembl_gene_id), ENSEMBL=as.character(query$ensembl_gene_id), GENES=as.character(query$external_gene_name), stringsAsFactors=FALSE)
    query <- cbind(query, paralog_perc_id=100, paralog_perc_id_r1=100)
    result <- result[nchar(result$hsapiens_paralog_ensembl_gene) > 0,] # removes rows with no paralog matches
    paralog <- data.frame(QueryID=result$ensembl_gene_id, ENSEMBL=result$hsapiens_paralog_ensembl_gene, GENES=result$hsapiens_paralog_associated_gene_name, paralog_perc_id=result$hsapiens_paralog_perc_id, paralog_perc_id_r1=result$hsapiens_paralog_perc_id_r1, stringsAsFactors=FALSE)
    resultDF <- rbind(query, paralog)
    resultDF <- resultDF[order(resultDF$QueryID),]

    ## Retrieve for ENSEMBL gene IDs in results the corresponding UniProt IDs.
    ## Note, this needs to be done in separate query since attributes from multiple
    ## Attribue pages are not allowed.
    uniprot <- getBM(attributes = c("ensembl_gene_id",                                                                                                                        
                                    "uniprot_gn_symbol",
                                    "uniprotswissprot",
                                    "uniprotsptrembl",
                                    "description"),
                    filters = "ensembl_gene_id",                                                                                                                          
                    values = resultDF$ENSEMBL,                                                                                                                                        
                    mart = mart)                

    ## For UniProtKB/Swiss-Prot IDs (curated entries)
    up_sp <- tapply(uniprot$uniprotswissprot, uniprot$ensembl_gene_id, function(x) unique(as.character(x)[nchar(x)>0]), simplify=FALSE)
    index <- sapply(up_sp, length)
    up_sp[index==0] <- ""
    index[index==0] <- 1
    index <- index[as.character(resultDF$ENSEMBL)]
    up_sp <- up_sp[names(index)] 
    resultDF <- cbind(resultDF[rep(seq_along(resultDF$ENSEMBL), index),], ID_up_sp=as.character(unlist(up_sp)))
    
    ## UniProtKB/TrEMBL IDs (less curated entries)
    up_sp_tr <- tapply(uniprot$uniprotsptrembl, uniprot$ensembl_gene_id, function(x) unique(as.character(x)[nchar(x)>0]), simplify=FALSE)
    index <- sapply(up_sp_tr, length)
    up_sp_tr[index==0] <- ""
    index[index==0] <- 1
    index <- index[as.character(resultDF$ENSEMBL)]
    up_sp_tr <- up_sp_tr[names(index)] 
    resultDF <- cbind(resultDF[rep(seq_along(resultDF$ENSEMBL), index),], ID_up_sp_tr=as.character(unlist(up_sp_tr)))
    row.names(resultDF) <- NULL

    ## Assemble result list
    IDMresult[IDMresult == ""] <- NA # This makes behavior more similar to result from getUniprotIDs
    for(i in colnames(IDMresult)) if(is.factor(IDMresult[,i])) IDMresult[,i] <- as.character(IDMresult[,i])
    resultDF[resultDF == ""] <- NA
    for(i in colnames(resultDF)) if(is.factor(resultDF[,i])) resultDF[,i] <- as.character(resultDF[,i])
    resultList <- list(IDM=IDMresult, SSNN=resultDF)
    return(resultList)
}
## Usage:
# queryBy <- list(molType="gene", idType="external_gene_name", ids=c("ZPBP", "MAPK1", "EGFR"))
# queryBy <- list(molType="gene", idType="ensembl_gene_id", ids=c("ENSG00000042813", "ENSG00000100030", "ENSG00000146648"))
# result <- getParalogs(queryBy)
# result$IDM[1:4,]
# result$SSNN[1:4,]
# sapply(result, dim, simplify=FALSE)
# sapply(names(result), function(x) unique(na.omit(result[[x]]$ID_up_sp)))

#####################################################
## Function to query known drug-target annotations ##
#####################################################

## Function to query known drug-target annotations
drugTargetAnnot <- function(queryBy=list(molType=NULL, idType=NULL, ids=NULL), cmpid_file=file.path(config$resultsPath,"cmp_ids.rds"), config=genConfig()) {
    if(any(names(queryBy) != c("molType", "idType", "ids"))) stop("All three list components in 'queryBy' (named: 'molType', 'idType' and 'ids') need to be present.")
    if(any(sapply(queryBy, length) == 0)) stop("All components in 'queryBy' list need to be populated with corresponding character vectors.")

    
    ## Load ChEMBL SQLite
    ## ChEMBL SQLite downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
    ## after unpacking you get chembl_xx.db
	 dbpath = config$chemblDbPath
    mydb <- dbConnect(SQLite(), dbpath)
    
    ## Input file for following step was generated by cmpIdMapping()
    cmp_ids <- readRDS(cmpid_file)
    rownames(cmp_ids) <- as.character(cmp_ids$molregno)
    
    ## Query bioactivity annotations
        ## Kevin's version
        # mainQuery <- "SELECT d.molregno, e.chembl_id, e.pref_name AS Drug_Name, d.mechanism_of_action AS MOA, d.action_type AS Action_Type, 
        #                                           e.first_approval AS First_Approval, a.chembl_id AS ChEMBL_TID, d.tid AS TID, c.accession AS UniProt_ID, 
        #                                           c.description AS Desc, c.organism AS Organism,
        #                                           GROUP_CONCAT(DISTINCT (f.mesh_id || ': ' || f.mesh_heading)) AS Mesh_Indication
        #                          FROM target_dictionary AS a
        #                                         JOIN target_components AS b USING(tid)
        #                                         JOIN component_sequences AS c USING(component_id)
        #                                         JOIN drug_mechanism AS d USING(tid)
        #                                         JOIN molecule_dictionary AS e USING(molregno)
        #                                         LEFT JOIN drug_indication AS f USING(molregno)"
    
    ## Query by targets
    if(queryBy$molType=="protein" & grepl("Uniprot", queryBy$idType[[1]], ignore.case=TRUE)) {
        idvec <- paste0("(\"", paste(queryBy$ids, collapse="\", \""), "\")") 
        
        ## Old version: delete
        # myquery <- dbSendQuery(mydb, paste("SELECT drug_mechanism.molregno, molecule_dictionary.chembl_id, 
        #                                             drug_mechanism.mechanism_of_action, drug_mechanism.tid, 
        #                                             target_components.component_id, component_sequences.accession AS UniProt_ID, 
        #                                             component_sequences.description, component_sequences.organism,
        #                                             drug_indication.mesh_id, drug_indication.mesh_heading
        #                           FROM drug_mechanism, molecule_dictionary, target_components, component_sequences, drug_indication
        #                           ON drug_mechanism.molregno = molecule_dictionary.molregno
        #                             AND drug_mechanism.tid = target_components.tid
        #                             AND target_components.component_id = component_sequences.component_id
        #                             AND molecule_dictionary.molregno = drug_indication.molregno
        #                             WHERE component_sequences.accession IN", idvec))
        
        ## Kevin's version
        # myquery <- dbSendQuery(mydb, paste(mainQuery,
        #                                     "WHERE c.accession IN", idvec,
        #                          "GROUP BY d.molregno, c.accession
        #                           ORDER BY c.accession, c.description"))
        
        ## Current version
        myquery <- dbSendQuery(mydb, paste("SELECT d.molregno, e.chembl_id, e.pref_name AS Drug_Name, d.mechanism_of_action AS MOA, d.action_type AS Action_Type, 
                                                   e.first_approval AS First_Approval, a.chembl_id AS ChEMBL_TID, d.tid AS TID, c.accession AS UniProt_ID, 
                                                   c.description AS Desc, c.organism AS Organism,
                                                   GROUP_CONCAT(DISTINCT (f.mesh_id || ': ' || f.mesh_heading)) AS Mesh_Indication
                                  FROM target_dictionary AS a
                                  LEFT JOIN target_components AS b ON a.tid = b.tid 
                                  LEFT JOIN component_sequences AS c ON b.component_id = c.component_id 
                                  LEFT JOIN drug_mechanism AS d ON a.tid = d.tid 
                                  LEFT JOIN molecule_dictionary AS e ON d.molregno = e.molregno 
                                  LEFT JOIN drug_indication AS f ON e.molregno = f.molregno
                                  WHERE c.accession IN", idvec,
                                  "GROUP BY d.molregno, c.accession
                                   ORDER BY c.accession, c.description"))
        activityassays <- dbFetch(myquery)
        dbClearResult(myquery)
    }
    
    ## Query by compounds
    if(queryBy$molType=="cmp") {
        ## ID conversions
        if(queryBy$idType=="molregno") {
            cmpvec <- queryBy$ids 
        }
        if(queryBy$idType=="chembl_id") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$chembl_id),]
            cmpvec <- as.character(cmp_ids$molregno); 
				names(cmpvec) <- as.character(cmp_ids$chembl_id)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        if(queryBy$idType=="PubChem_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$PubChem_ID),]
            cmpvec <- as.character(cmp_ids$molregno); 
				names(cmpvec) <- as.character(cmp_ids$PubChem_ID)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        if(queryBy$idType=="DrugBank_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$DrugBank_ID),]
            cmpvec <- as.character(cmp_ids$molregno); 
				names(cmpvec) <- as.character(cmp_ids$DrugBank_ID)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        idvec <- paste0("(\"", paste(cmpvec, collapse="\", \""), "\")") 
        
        ## Old version: delete
        # myquery <- dbSendQuery(mydb, paste("SELECT drug_mechanism.molregno, molecule_dictionary.chembl_id, molecule_dictionary.pref_name, 
        #                                             drug_mechanism.mechanism_of_action, drug_mechanism.action_type, drug_mechanism.tid, 
        #                                             target_components.component_id, target_dictionary.chembl_id AS ChEMBL_TID, component_sequences.accession AS UniProt_ID, 
        #                                             component_sequences.description AS Desc, component_sequences.organism AS Organism,
        #                                             drug_indication.mesh_id, drug_indication.mesh_heading
        #                           FROM drug_mechanism, molecule_dictionary, target_dictionary, target_components, component_sequences, drug_indication
        #                           ON drug_mechanism.molregno = molecule_dictionary.molregno
        #                             AND drug_mechanism.tid = target_components.tid
        #                             AND drug_mechanism.tid = target_dictionary.tid
        #                             AND target_components.component_id = component_sequences.component_id
        #                             AND molecule_dictionary.molregno = drug_indication.molregno
        #                             WHERE molecule_dictionary.molregno IN", idvec))
        
        ## Kevin's version
        # myquery <- dbSendQuery(mydb, paste(mainQuery,
        #                          "WHERE e.molregno IN", idvec,
        #                          "GROUP BY e.molregno, c.accession
        #                           ORDER BY e.molregno, e.pref_name"))
        

        ## Current version
        myquery <- dbSendQuery(mydb, paste("SELECT a.molregno, a.chembl_id, a.pref_name AS Drug_Name, b.mechanism_of_action AS MOA, 
																	b.action_type AS Action_Type, 
                                                   a.first_approval AS First_Approval,       
                                                   c.chembl_id AS ChEMBL_TID, b.tid AS TID, e.accession AS UniProt_ID, 
                                                   e.description AS Desc, e.organism AS Organism,
                                                   GROUP_CONCAT(DISTINCT (f.mesh_id || ': ' || f.mesh_heading)) AS Mesh_Indication
                                  FROM molecule_dictionary AS a
                                  LEFT JOIN drug_mechanism AS b ON a.molregno = b.molregno 
                                  LEFT JOIN target_dictionary AS c ON b.tid = c.tid
                                  LEFT JOIN target_components AS d ON b.tid = d.tid
                                  LEFT JOIN component_sequences AS e ON d.component_id = e.component_id
                                  LEFT JOIN drug_indication AS f ON a.molregno = f.molregno
                                  WHERE a.molregno IN", idvec,
                                  "GROUP BY a.molregno, e.accession
                                   ORDER BY a.molregno, a.pref_name"))
        activityassays <- dbFetch(myquery)
        dbClearResult(myquery)
    }
    resultDF <- data.frame(cmp_ids[as.character(activityassays$molregno),-c(2,4)], activityassays[,-c(1,2)])
    dbDisconnect(mydb)
    
    ## Remove rows with identical values in all fields
    resultDF <- resultDF[!duplicated(apply(resultDF, 1, paste, collapse="_")),]
    ## Add query column and sort rows in result table according to query
    index_list <- sapply(queryBy$ids, function(x) which(toupper(resultDF[,queryBy$idType]) %in% toupper(x)), simplify=FALSE)
    index_list[sapply(index_list,length)==0] <- Inf
    index_df <- data.frame(ids=rep(names(index_list), sapply(index_list, length)), rowids=unlist(index_list))
    resultDF <- data.frame(QueryIDs=index_df[,1], resultDF[as.numeric(index_df$rowids),])
    rownames(resultDF) <- NULL
    return(resultDF)
}

#########################################
## Query Known Drug-Target Annotations ##
#########################################
## Function to generate drugTargetAnnot.xls table
drugTargetAnnotTable <- function(outfile, rerun=TRUE,config=genConfig()) {
    if(rerun==TRUE) {
        ## Load ChEMBL SQLite
        ## ChEMBL SQLite downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
        ## after unpacking you get chembl_24.db
        chembldb <- config$chemblDbPath
        mydb <- dbConnect(SQLite(), chembldb)

        ## Joins for drug mech
        allmech <- dbGetQuery(mydb, 'SELECT * FROM drug_mechanism, molecule_dictionary WHERE drug_mechanism.molregno = molecule_dictionary.molregno')
        
        ## Set values in "mechanism_comment" column to NA since it contains weird characters resulting in row duplications when writing to files
        allmech[ , colnames(allmech) %in% "mechanism_comment"] <- rep(NA, nrow(allmech))

        ## Get UniChem Info
        ## ChEMBL to DrubBank mapping in src1src2.txt.gz: 
        chembl2drugbank <- read.delim(gzfile(file.path(config$downloadPath,"src1src2.txt.gz")))
        chembl2drugbank_vec <- tapply(as.character(chembl2drugbank[,2]), factor(chembl2drugbank[,1]), paste, collapse=", ")
        ## ChEMBL to PubChem CID mapping in src1src22.txt.gz
        chembl2pubchem <- read.delim(gzfile(file.path(config$downloadPath,"src1src22.txt.gz")))
        chembl2pubchem_vec <- tapply(as.character(chembl2pubchem[,2]), factor(chembl2pubchem[,1]), paste, collapse=", ")
        ## ChEMBL to ChEBI mapping in src1src7.txt.gz
        chembl2chebi <- read.delim(gzfile(file.path(config$downloadPath,"src1src7.txt.gz")))
        chembl2chebi_vec <- tapply(as.character(chembl2chebi[,2]), factor(chembl2chebi[,1]), paste, collapse=", ")
        allmech <- cbind(allmech, PubChem_ID=chembl2pubchem_vec[allmech$chembl_id], DrugBank_ID=chembl2drugbank_vec[allmech$chembl_id], ChEBI_ID=chembl2chebi_vec[allmech$chembl_id])

        ## Joins for target ids
        targetids <- dbGetQuery(mydb, 'SELECT * FROM target_components, component_sequences WHERE target_components.component_id = component_sequences.component_id')
        tid2accession <- tapply(targetids$accession, as.factor(targetids$tid), paste, collapse=", ")
        tid2description <- tapply(targetids$description, as.factor(targetids$tid), paste, collapse=", ")
        tid2organism <- tapply(targetids$organism, as.factor(targetids$tid), paste, collapse=", ")
        target_dict <- dbGetQuery(mydb, 'SELECT * FROM target_dictionary')
        chembl_tid <- as.character(target_dict$chembl_id); names(chembl_tid) <- as.character(target_dict$tid)
        tidDF <- data.frame(row.names=names(tid2accession), ChEMBL_TID=chembl_tid[names(tid2accession)], UniProt_ID=tid2accession, Desc=tid2description[names(tid2accession)], Organism=tid2organism[names(tid2accession)])
        
        ## Join mesh and efo indication
        indication <- dbGetQuery(mydb, 'SELECT * FROM drug_indication')
        mesh <- paste(indication$mesh_id, indication$mesh_heading, sep=": ")
        mesh2molregno <- tapply(mesh, as.factor(indication$molregno), paste, collapse=", ")
        efo <- paste(indication$efo_id, indication$efo_term, sep=": ")
        efo2molregno <- tapply(efo, as.factor(indication$molregno), paste, collapse=", ")

        ## Assemble everything
        allmech <- data.frame(row.names=seq_along(allmech[,1]), 
                              allmech, 
                              tidDF[as.character(allmech$tid),],
                              MeshIndication=mesh2molregno[as.character(allmech$molregno)],
                              EFOIndication=efo2molregno[as.character(allmech$molregno)])
        write.table(allmech, outfile, row.names=FALSE, quote=FALSE, sep="\t")
    }
}
## Usage:
# drugTargetAnnotTable(outfile="results/drugTargetAnnot.xls", rerun=FALSE)

## Compare with online drug-target table from ChEMBL obtained from manual download here: https://www.ebi.ac.uk/chembl/old/drug/targets

## Query functions for drugTargetAnnot.xls
getDrugTarget <- function(dt_file=file.path(config$resultsPath,"drugTargetAnnot.xls"), queryBy=list(molType=NULL, idType=NULL, ids=NULL), 
								  id_mapping=c(chembl="chembl_id", pubchem="PubChem_ID", uniprot="UniProt_ID"), columns,config=genConfig()) {
	#stop("dt_file: ",dt_file)
    dt_file <- read.delim(dt_file)
    myid <- id_mapping[queryBy$idType]
    .queryFct <- function(dt_file, myid) {
        idlist <- strsplit(as.character(dt_file[,myid]), ", ")
        index_list <- sapply(queryBy$ids, function(x) which(sapply(seq_along(idlist), 
                                          function(y) any(tolower(idlist[[y]]) %in% tolower(x)))), simplify=FALSE)
        ## To include no matches in result, inject Inf in corresponding slots 
        index_list[sapply(index_list,length)==0] <- Inf
        index_df <- data.frame(ids=rep(names(index_list), sapply(index_list, length)), rowids=unlist(index_list))
        df <- data.frame(QueryIDs=index_df[,1], dt_file[as.numeric(index_df$rowids),])
        return(df)
    }
    df <- .queryFct(dt_file, myid)
    dfsub <- df[,columns]
    rownames(dfsub) <- NULL
    return(dfsub)
}
## Usage:
# id_mapping <- c(chembl="chembl_id", pubchem="PubChem_ID", uniprot="UniProt_ID", drugbank="DrugBank_ID")
# queryBy <- list(molType="cmp", idType="chembl", ids=c("CHEMBL25", "CHEMBL1742471"))
# getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52)) 
# queryBy <- list(molType="cmp", idType="pubchem", ids=c("2244", "65869", "2244"))
# getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52)) 
# queryBy <- list(molType="protein", idType="uniprot", ids=c("P43166", "P00915", "P43166"))
# getDrugTarget(dt_file="results/drugTargetAnnot.xls", queryBy=queryBy, id_mapping, columns=c(1,5,8,16,17,39,46:52))

#############
## DrugAge ##
#############
## Source file of DrugAge is linked here: http://genomics.senescence.info/drugs/ 
processDrugage <- function(drugagefile=file.path(config$resultsPath,"drugage_id_mapping.xls"), redownloaddrugage=TRUE,config=genConfig()) {
    if(redownloaddrugage==TRUE) {
        download.file("http://genomics.senescence.info/drugs/dataset.zip", file.path(config$downloadPath,"dataset.zip"))
        unzip(file.path(config$downloadPath,"dataset.zip"), exdir=config$downloadPath)
    }
    drugage <- read.csv(file.path(config$downloadPath,"drugage.csv"))
    drugage <- drugage[!duplicated(drugage$compound_name),]
    ## Get mol_dict table from chembl_db
    chembldb <- config$chemblDbPath
    mydb <- dbConnect(SQLite(), chembldb)
    mol_dict <- dbGetQuery(mydb, 'SELECT pref_name,chembl_id FROM molecule_dictionary')
    mol_dict_sub <- mol_dict[!is.na(mol_dict$pref_name),] # mol_dict from below
    prefname <- as.character(mol_dict_sub$pref_name)
    chemblid <- as.character(mol_dict_sub$chembl_id)
    fact <- tapply(chemblid, factor(prefname), paste, collapse=", ")
    ## Load ChEMBL to PubChem CID and DrugBank ID mappings (generated above with downloadUniChem Fct)
    chembl2pubchem <- read.delim(gzfile(file.path(config$downloadPath,"src1src22.txt.gz")))
    chembl2pubchem_vec <- tapply(as.character(chembl2pubchem[,2]), factor(chembl2pubchem[,1]), paste, collapse=", ")
    chembl2drugbank <- read.delim(gzfile(file.path(config$downloadPath,"src1src2.txt.gz")))
    chembl2drugbank_vec <- tapply(as.character(chembl2drugbank[,2]), factor(chembl2drugbank[,1]), paste, collapse=", ")
    ## Assemble results
    drugage <- cbind(drugage, pref_name=names(fact[toupper(drugage$compound_name)]), chembl_id=fact[toupper(drugage$compound_name)])
    drugage <- cbind(drugage, pubchem_cid=chembl2pubchem_vec[as.character(drugage$chembl_id)], drugbank_id=chembl2drugbank_vec[as.character(drugage$chembl_id)])
    write.table(drugage, drugagefile, row.names=FALSE, quote=FALSE, sep="\t")
}
## Usage:
# processDrugage(drugagefile="results/drugage_id_mapping.xls", redownloaddrugage=FALSE)
## Now the missing IDs need to be added manually. A semi-manual approach is to use this web service: https://cts.fiehnlab.ucdavis.edu/batch

########################################################
## Integration with Therapeutic Target Database (TTD) ##
########################################################
## Query performed for Nik Schork
## ID crossmatching table
transformTTD <- function(ttdfile=file.path(config$downloadPath,"TTD_IDs.txt"), redownloadTTD=TRUE,config=genConfig()) {
    if(redownloadTTD==TRUE) download.file("https://db.idrblab.org/ttd/sites/default/files/ttd_database/P1-02-TTD_crossmatching.txt", ttdfile)
    df <- read.delim(ttdfile, skip=12, header=FALSE, col.names=c("TTD_ID", "No", "Source", "Name"))
    df_sub <- df[grep("DrugName|PubChem CID", df$Source),]
    list_sub <- split(df_sub, df_sub$TTD_ID)
    list_sub <- list_sub[lapply(list_sub, nrow)==2]
    list_sub <- lapply(list_sub, t)
    list_sub <- lapply(list_sub, tail, 1)
    df_sub <- do.call("rbind", list_sub)
    row.names(df_sub) <- names(list_sub)
    df_sub <- data.frame(TTD_ID=rownames(df_sub), CMP_Name=df_sub[,1], PubChem_ID=gsub("^CID ", "", as.character(df_sub[,2])))
    return(df_sub)
}
## Usage:
# df_sub <- transformTTD(ttdfile="downloads/TTD_IDs.txt", redownloadTTD=FALSE) 
# df_sub[1:4,]


############################
## Query Bioactivity Data ##
############################
## Function to generate CMP ID mappings UniChem (current version very slow)
cmpIdMapping <- function(outfile=file.path(config$resultsPath,"cmp_ids.rds"), rerun=TRUE,config=genConfig()) {
    if(rerun==TRUE) {
        ## Load ChEMBL SQLite
        ## ChEMBL SQLite downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
        ## after unpacking you get chembl_24.db
        chembldb <- config$chemblDbPath
        mydb <- dbConnect(SQLite(), chembldb)

        ## Drug ID mappings
        chemblid_lookup <- dbGetQuery(mydb, 'SELECT * FROM chembl_id_lookup WHERE "entity_type" = "COMPOUND"')
        colnames(chemblid_lookup)[colnames(chemblid_lookup)=="entity_id"] <- "molregno"
       
        ## CMP ID mappings from UniChem (very slow!)
        ## ChEMBL to DrubBank mapping in src1src2.txt.gz: 
        chembl2drugbank <- read.delim(gzfile(file.path(config$downloadPath,"src1src2.txt.gz")))
        chembl2drugbank_vec <- tapply(as.character(chembl2drugbank[,2]), factor(chembl2drugbank[,1]), paste, collapse=", ")
        ## ChEMBL to PubChem CID mapping in src1src22.txt.gz
        chembl2pubchem <- read.delim(gzfile(file.path(config$downloadPath,"src1src22.txt.gz")))
        chembl2pubchem_vec <- tapply(as.character(chembl2pubchem[,2]), factor(chembl2pubchem[,1]), paste, collapse=", ")
        ## ChEMBL to ChEBI mapping in src1src7.txt.gz
        chembl2chebi <- read.delim(gzfile(file.path(config$downloadPath,"src1src7.txt.gz")))
        chembl2chebi_vec <- tapply(as.character(chembl2chebi[,2]), factor(chembl2chebi[,1]), paste, collapse=", ")
        ## Create final CMP ID mapping table (very slow!!)
        # cmp_ids <- cbind(chemblid_lookup, PubChem_ID=chembl2pubchem_vec[as.character(chemblid_lookup$chembl_id)], DrugBank_ID=chembl2drugbank_vec[as.character(chemblid_lookup$chembl_id)], ChEBI_ID=chembl2chebi_vec[as.character(chemblid_lookup$chembl_id)])
        cmp_ids <- cbind(chemblid_lookup, 
								 PubChem_ID=chembl2pubchem_vec[as.character(chemblid_lookup$chembl_id)], 
								 DrugBank_ID=chembl2drugbank_vec[as.character(chemblid_lookup$chembl_id)])
        saveRDS(cmp_ids, outfile)
    }
}
## Usage: 
# cmpIdMapping(outfile="results/cmp_ids.rds", rerun=FALSE)
# cmp_ids <- readRDS("results/cmp_ids.rds")

## Function to query bioactivity data by target or compound ids
drugTargetBioactivity <- function( queryBy=list(molType=NULL, idType=NULL, ids=NULL), cmpid_file=file.path(config$resultsPath,"cmp_ids.rds"),config=genConfig()) {
    if(any(sapply(queryBy, is.null))) stop("All components in 'queryBy' list need to be populated with corresponding character vectors.")
    
    ## Load ChEMBL SQLite
    ## ChEMBL SQLite downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
    ## after unpacking you get chembl_xx.db
	 dbpath = config$chemblDbPath
    mydb <- dbConnect(SQLite(), dbpath)
    
    ## Input file for following step was generated by cmpIdMapping()
    cmp_ids <- readRDS(cmpid_file)
    rownames(cmp_ids) <- as.character(cmp_ids$molregno)
    
    ## Query bioactivity
    
    ## Kevin's version
    # mainQuery <- "SELECT activities.molregno, pref_name, activity_id, assays.chembl_id AS chembl_assay_id, accession AS UniProt_ID, 
    #                             component_sequences.description, organism, activities.standard_value, activities.standard_units, 
    #                             activities.standard_flag, activities.standard_type
    #                    FROM molecule_dictionary
    #                          JOIN activities USING(molregno)
    #                          JOIN assays USING(assay_id)
    #                          JOIN target_components USING(tid)
    #                          JOIN component_sequences USING(component_id)"
    
    ## Query by targets
    if(queryBy$molType=="protein" & queryBy$idType=="uniprot") {
        idvec <- paste0("(\"", paste(queryBy$ids, collapse="\", \""), "\")") 
        
        ## Kevin's version
        # myquery <- dbSendQuery(mydb, paste(mainQuery, "WHERE component_sequences.accession IN", idvec))
        

        ## Current version
        #myquery <- dbSendQuery(mydb, paste("SELECT activities.molregno, pref_name, activity_id, assays.chembl_id AS chembl_assay_id, accession AS UniProt_ID, component_sequences.description, organism, drug_indication.mesh_heading, activities.standard_value, activities.standard_units, activities.standard_flag, activities.standard_type
        myquery <- dbSendQuery(mydb, paste("SELECT activities.molregno, pref_name, activity_id, assays.chembl_id AS chembl_assay_id, accession AS UniProt_ID, component_sequences.description, organism, activities.standard_value, activities.standard_units, activities.standard_flag, activities.standard_type
                                  FROM activities, assays, target_components, component_sequences, molecule_dictionary
                                  ON activities.assay_id = assays.assay_id
                                    AND assays.tid = target_components.tid
                                    AND target_components.component_id = component_sequences.component_id
                                    AND activities.molregno = molecule_dictionary.molregno
                                    WHERE component_sequences.accession IN", idvec))
        activityassays <- dbFetch(myquery)
        dbClearResult(myquery)
    }
    
    ## Query bioactivity by compounds
    if(queryBy$molType=="cmp") {
        ## ID conversions
        if(queryBy$idType=="molregno") {
            cmpvec <- queryBy$ids 
        }
        if(queryBy$idType=="chembl_id") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$chembl_id),]
            cmpvec <- as.character(cmp_ids$molregno); names(cmpvec) <- as.character(cmp_ids$chembl_id)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        if(queryBy$idType=="PubChem_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$PubChem_ID),]
            cmpvec <- as.character(cmp_ids$molregno); names(cmpvec) <- as.character(cmp_ids$PubChem_ID)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        if(queryBy$idType=="DrugBank_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$DrugBank_ID),]
            cmpvec <- as.character(cmp_ids$molregno); names(cmpvec) <- as.character(cmp_ids$DrugBank_ID)     
            cmpvec <- cmpvec[queryBy$ids]
        }
        idvec <- paste0("(\"", paste(cmpvec, collapse="\", \""), "\")") 
        
        ## Kevin's version
        # myquery <- dbSendQuery(mydb, paste(mainQuery, "WHERE activities.molregno IN", idvec))

        ## Current version
        myquery <- dbSendQuery(mydb, paste("SELECT activities.molregno, pref_name, activity_id, assays.chembl_id AS chembl_assay_id, accession AS UniProt_ID, component_sequences.description, organism, activities.standard_value, activities.standard_units, activities.standard_flag, activities.standard_type
                                  FROM activities, assays, target_components, component_sequences, molecule_dictionary
                                  ON activities.assay_id = assays.assay_id
                                    AND assays.tid = target_components.tid
                                    AND target_components.component_id = component_sequences.component_id
                                    AND activities.molregno = molecule_dictionary.molregno
                                    WHERE activities.molregno IN", idvec))
        activityassays <- dbFetch(myquery)
        dbClearResult(myquery)
    }
    resultDF <- data.frame(cmp_ids[as.character(activityassays$molregno),-c(2,4)], activityassays[,-1])
    dbDisconnect(mydb)
    rownames(resultDF) <- NULL
    return(resultDF)
}

## Usage:
# chembldb <- "/bigdata/girkelab/shared/lcshared/chemoinformatics/compoundDBs/chembl_24/chembl_24_sqlite/chembl_24.db"
# queryBy <- list(molType="protein", idType="uniprot", ids=c("P05979", "P35354", "P33033", "Q8VCT3", "P29475", "P51511"))
# qresult <- drugTargetBioactivity( queryBy, cmpid_file="results/cmp_ids.rds")
# qresult[1:4,]
# queryBy <- list(molType="cmp", idType="molregno", ids=c("101036", "101137", "1384464")) 
# qresult <- drugTargetBioactivity( queryBy, cmpid_file="results/cmp_ids.rds")
# qresult[1:4,]
# queryBy <- list(molType="cmp", idType="DrugBank_ID", ids=c("DB00945", "DB00316", "DB01050")) 
# qresult <- drugTargetBioactivity(queryBy, cmpid_file="results/cmp_ids.rds")
# qresult[1:4,]
# queryBy <- list(molType="cmp", idType="PubChem_ID", ids=c("2244", "3672", "1983")) 
# qresult <- drugTargetBioactivity( queryBy, cmpid_file="results/cmp_ids.rds")
# qresult[1:4,]

#################################
## Gene to Protein ID Mappings ##
#################################
## Return for gene or protein IDs a mapping table containing: ENSEMBL Gene IDs,
## Gene Names/Symbols, UniProt IDs and ENSEMBL Protein IDs. The results are returned
## in a list where the first slot contains the ID mapping table and the subsequent slots
## the corresponding named character vectors for: ens_gene_id, up_ens_id, and up_gene_id.
## Currently, the query IDs can be Gene Names, Ensembl Gene IDs and UniProt IDs.
getSymEnsUp <- function(EnsDb="EnsDb.Hsapiens.v86", ids, idtype) {
    
	 require(EnsDb,character.only=TRUE)
    ## ID columns to return
    idcolumns <- c("gene_id", "gene_name", "uniprot_id", "protein_id") 

    ## With gene names
    if(idtype=="GENE_NAME") {
        idDF <- genes(get(EnsDb), filter = GeneNameFilter(ids), columns = idcolumns) 
        idDF <- idDF[!is.na(values(idDF)$uniprot_id),]
        idDF <- values(idDF)
    ## With ENSEBML Gene IDs
    } else if(idtype=="ENSEMBL_GENE_ID") {
        idDF <- genes(get(EnsDb), filter = GeneIdFilter(ids), columns = idcolumns) 
        idDF <- idDF[!is.na(values(idDF)$uniprot_id),]
        idDF <- values(idDF)
    ## With UniProt IDs
    } else if(idtype=="UNIPROT_ID") {
        idDF <- proteins(get(EnsDb), filter = UniprotFilter(ids), columns = idcolumns)
        idDF <- idDF[!is.na(idDF$uniprot_id),]
    } else {
        stop("Argument 'idtype' can only be assigned one of: 'GENE_NAME', 'ENSEMBL_GENE_ID', 'UNIPROT_ID'.")
    }
    ## Assemble ID mapping table
    idDF <- idDF[,idcolumns]
    row.names(idDF) <- NULL

    ## Assemble named char vector: ens_gene_id
    idDF <- as.data.frame(idDF)
    ens_gene_idDF <- idDF[!duplicated(idDF$gene_id),]
    ens_gene_id <- as.character(ens_gene_idDF$gene_name)
    names(ens_gene_id) <- as.character(ens_gene_idDF$gene_id)
    
    ## Assemble named char vector: up_ens_id
    up_ens_idDF <- idDF[!duplicated(idDF$uniprot_id),]
    up_ens_id <- as.character(up_ens_idDF$gene_id)
    names(up_ens_id) <- as.character(up_ens_idDF$uniprot_id)
    
    ## Assemble named char vector: up_gene_id
    up_gene_idDF <- idDF[!duplicated(idDF$uniprot_id),]
    up_gene_id <- as.character(up_gene_idDF$gene_name)
    names(up_gene_id) <- as.character(up_gene_idDF$uniprot_id)
   
    ## Return result as list
    returnList <- list(idDF=idDF, ens_gene_id=ens_gene_id, up_ens_id=up_ens_id, up_gene_id=up_gene_id)
    return(returnList)
}
## Usage:
# gene_name <- c("CA7", "CFTR")
# getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=gene_name, idtype="GENE_NAME")
# ensembl_gene_id <- c("ENSG00000001626", "ENSG00000168748")
# getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=ensembl_gene_id, idtype="ENSEMBL_GENE_ID")
# uniprot_id <- c("P43166", "P13569") 
# getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=uniprot_id, idtype="UNIPROT_ID")

##############################################
## Convenience Drug-Target Wrapper Function ##
##############################################
## Meta-function to obtain both drug-target annotation and bioassay data 
runDrugTarget_Annot_Bioassay <- function(res_list, up_col_id="ID", ens_gene_id,  cmpid_file=file.path(config$resultsPath,"cmp_ids.rds") ,config=genConfig(), ...) {
	## (1.1) Check input validity 
	## To be added, once class system has been introduced...
	## (1.2) Input ID mappings, note: Ensembl query IDs with similar sequences will often map to the same related/third Ensemble gene resulting in paralog genes/proteins mapping to several query genes. Thus, this one-to-many relationship needs to be resolved for both geneids and ensids.
	id_list <- sapply(names(res_list), function(x) unique(na.omit(as.character(res_list[[x]][,up_col_id]))), simplify=FALSE)
	ensids <- tapply(res_list[[2]]$QueryID, res_list[[2]][,up_col_id], function(x) as.character(unique(x)), simplify=FALSE) 
    geneids <- sapply(names(ensids), function(x) ens_gene_id[ensids[[x]]], simplify=FALSE)
    ensids <- sapply(names(ensids), function(x) paste(ensids[[x]], collapse=", "))
    geneids <- sapply(names(geneids), function(x) paste(geneids[[x]], collapse=", "))
    ensids <- unlist(ensids[nchar(names(ensids)) > 0]) 
    geneids <- unlist(geneids[nchar(names(geneids)) > 0]) 

	## (1.3) Obtain drug-target annotations for IDM
	idm_ids <- id_list$IDM
	if(length(idm_ids) > 0) {
	    queryBy <- list(molType="protein", idType="UniProt_ID", ids=idm_ids)
	    qresultIDM <- drugTargetAnnot( queryBy, cmpid_file=cmpid_file,config=config)
	    qresultIDM <- data.frame(GeneName=as.character(geneids[as.character(qresultIDM$QueryIDs)]),
	                             Ensembl_IDs=as.character(ensids[as.character(qresultIDM$QueryIDs)]), 
	                             qresultIDM)
    } else {
        qresultIDM <- data.frame()
    }

	## (1.4) Obtain drug-target annotations for SSNN excluding IDM overlaps (SSNN_noIDM)
	ssnn_ids <- id_list$SSNN[!id_list$SSNN %in% idm_ids] # UniProt IDs only in SSNN
	if(length(ssnn_ids) > 0) {
	    queryBy <- list(molType="protein", idType="UniProt_ID", ids=ssnn_ids)
	    qresultSSNN <- drugTargetAnnot( queryBy, cmpid_file=cmpid_file,config=config)
	    qresultSSNN <- data.frame(GeneName=as.character(geneids[as.character(qresultSSNN$QueryIDs)]), 
	                              Ensembl_IDs=as.character(ensids[as.character(qresultSSNN$QueryIDs)]), 
	                              qresultSSNN)
        ## Prepend "Query_" to indicate that returned genes IDs are often not the ones encoding SSNN UniProt IDs
        if(length(qresultSSNN[, "Ensembl_IDs"]) > 0) qresultSSNN[, "Ensembl_IDs"] <- paste0("Query_", qresultSSNN[, "Ensembl_IDs"])
        qresultSSNN[qresultSSNN[,"Ensembl_IDs"]=="Query_NA","Ensembl_IDs"] <- NA # Fixes Query_NA strings in NA cases
        if(length(qresultSSNN[, "GeneName"]) > 0 ) qresultSSNN[, "GeneName"] <- paste0("Query_", qresultSSNN[, "GeneName"])
        qresultSSNN[qresultSSNN[,"GeneName"]=="Query_NA","GeneName"] <- NA # Fixes Query_NA strings in NA cases
    } else {
        qresultSSNN <- data.frame()
    }

	## (1.5) Combine IDM and SSNN_noIDM results
	qresult <- rbind(cbind(ID_Mapping_Type=rep("IDM", nrow(qresultIDM)), qresultIDM), cbind(ID_Mapping_Type=rep("SSNN_noIDM", nrow(qresultSSNN)), qresultSSNN))
	qresult <- qresult[order(gsub("Query_", "", qresult$GeneName), qresult$ID_Mapping_Type),]
	colnames(qresult) <- c("ID_Mapping_Type", "GeneName", "Ensembl_IDs", "UniProt_QueryIDs", "CHEMBL_CMP_ID", "Molregno", "PubChem_ID", "DrugBank_ID", "Drug_Name", "MOA", "Action_Type", "First_Approval", "ChEMBL_TID", "TID", "UniProt_ID", "Target_Desc", "Organism", "Mesh_Indication")
	rownames(qresult) <- NULL

	## (2.1) Obtain bioassay data for IDM
	queryBy <- list(molType="protein", idType="uniprot", ids=idm_ids)
	qresultBAIDM <- drugTargetBioactivity( queryBy, cmpid_file=cmpid_file,config=config)
	qresultBAIDM <- data.frame(GeneName=as.character(geneids[as.character(qresultBAIDM$UniProt_ID)]), 
	                           Ensembl_IDs=as.character(ensids[as.character(qresultBAIDM$UniProt_ID)]), 
	                           qresultBAIDM)

	## (2.2) Obtain bioassay data for SSNN excluding IDM overlaps (SSNN_noIDM)
	queryBy <- list(molType="protein", idType="uniprot", ids=ssnn_ids)
	qresultBASSNN <- suppressWarnings(drugTargetBioactivity(queryBy, cmpid_file=cmpid_file,config=config))
	qresultBASSNN <- data.frame(GeneName=geneids[as.character(qresultBASSNN$UniProt_ID)], 
	                            Ensembl_IDs=ensids[as.character(qresultBASSNN$UniProt_ID)], 
	                            qresultBASSNN)
    ## Prepend "Query_" to indicate that returned genes IDs are often not the ones encoding SSNN UniProt IDs
    if(length(qresultBASSNN[, "Ensembl_IDs"]) > 0) qresultBASSNN[, "Ensembl_IDs"] <- paste0("Query_", qresultBASSNN[, "Ensembl_IDs"])
    qresultBASSNN[qresultBASSNN[,"Ensembl_IDs"]=="Query_NA","Ensembl_IDs"] <- NA # Fixes Query_NA strings in NA cases
    if(length(qresultBASSNN[, "GeneName"]) > 0) qresultBASSNN[, "GeneName"] <- paste0("Query_", qresultBASSNN[, "GeneName"])
    qresultBASSNN[qresultBASSNN[,"GeneName"]=="Query_NA","GeneName"] <- NA # Fixes Query_NA strings in NA cases

	## (2.3) Combine IDM and SSNN_noIDM results
	qresultBA <- rbind(cbind(ID_Mapping_Type=rep("IDM", nrow(qresultBAIDM)), qresultBAIDM), cbind(ID_Mapping_Type=rep("SSNN_noIDM", nrow(qresultBASSNN)), qresultBASSNN))
	qresultBA <- qresultBA[order(gsub("Query_", "", qresultBA$GeneName), qresultBA$ID_Mapping_Type),]
	colnames(qresultBA) <- c("ID_Mapping_Type", "GeneName", "Ensembl_IDs", "CHEMBL_CMP_ID", "Molregno", "PubChem_ID", "DrugBank_ID", "Drug_Name", "Activity_ID", "CHEMBL_Assay_ID", "UniProt_ID", "Target_Desc", "Organism", "Standard_Value", "Standard_Units", "Standard_Flag", "Standard_Type")
	rownames(qresultBA) <- NULL

	## (3.1) Export annotation and bioassay results in list
	drug_target_list <- list(Annotation=qresult, Bioassay=qresultBA)
	return(drug_target_list)
}
## Usage:
# ensembl_gene_id <- c("ENSG00000001626", "ENSG00000168748")
# idMap <- getSymEnsUp(EnsDb="EnsDb.Hsapiens.v86", ids=ensembl_gene_id, idtype="ENSEMBL_GENE_ID")
# ens_gene_id <- idMap$ens_gene_id
# res_list90 <- getUniprotIDs(taxId=9606, kt="ENSEMBL", keys=names(ens_gene_id), seq_cluster="UNIREF90", chunksize=1)
# queryBy <- list(molType="gene", idType="ensembl_gene_id", ids=names(ens_gene_id))
# res_list <- getParalogs(queryBy)
# drug_target_list <- runDrugTarget_Annot_Bioassay(res_list=res_list, up_col_id="ID_up_sp", ens_gene_id, chembldb, cmp_ids.rds) 



