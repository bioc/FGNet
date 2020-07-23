
######################################################
# refList
# buildGeneSets
# findGenesInTerm - TODO - 

######################################################
# geneIDtype <- "ENSEMBL"
# dbPackage <- "hgu133plus2.db" # organism or chip database
# if(is.null(genesUniverse)) genesUniverse <- refList(organism, geneIDtype)
refList <- function(dbPackage, geneIdType)
{
    if(!loadInstPkg(dbPackage)) stop(paste("Package", dbPackage, "is not available."))

    pkg.db <- eval(parse(text=dbPackage))
    AnnotationDbi::columns(pkg.db)
    
    if(!geneIdType %in% AnnotationDbi::columns(pkg.db)) stop(paste("geneIdType not available for ",dbPackage,". \nAvailable columns: ",paste(columns(pkg.db), collapse=", "), sep=""))
    
    return(AnnotationDbi::keys(pkg.db, keytype=geneIdType))
}


# organismDb <- "org.Hs.eg.db"
# geneIDtype <- "ENSEMBL"

# Builds gene-term sets based on Bioconductor databases: Organisms and GO/KEGG/REACTOME
# used by query_topGO and query_gage
buildGeneSets <- function(organismDb, geneIDtype, annotations=c("GO_BP","GO_MF","GO_CC","KEGG","REACTOME"), evidence=NULL, termLabel=c("TERM","ID"))
{
    if(termLabel[1] == "TERM") termAsLabel <- TRUE
    # Check and load required libraries
    if("REACTOME" %in% annotations) 
    {
        if(!loadInstPkg("reactome.db")) 
        {
            warning("Package 'reactome.db' is required to build REACTOME sets.")
            annotations <- annotations[which(!annotations%in%"REACTOME")]
        }
    }

    if(!loadInstPkg(organismDb)) stop("The organism package is not available.")

    org.db <- eval(parse(text=organismDb))
    
    if(!geneIDtype %in% AnnotationDbi::columns(org.db)) stop(paste("geneIDtype not available for ",organismDb,". \nAvailable columns: ",paste(columns(org.db), collapse=", "), sep=""))

    allGenes <- AnnotationDbi::keys(org.db, keytype=geneIDtype)
    allGenes <- gsub(pattern="\"", "",allGenes)
    ret <- list(genes2terms=NULL, terms2genes=NULL, asTable=NULL)
    #### Build GO tables:
    if(any(c("GO_BP","GO_MF","GO_CC") %in% annotations))
    {
        tableGoGenes <- suppressWarnings(AnnotationDbi::select(org.db, keys=allGenes, columns="GO", keytype=geneIDtype))
        tableGoGenes <- tableGoGenes[!is.na(tableGoGenes[,"GO"]),] # Bueno o malo?
        if(!is.null(evidence)) tableGoGenes <- tableGoGenes[tableGoGenes$EVIDENCE %in% evidence, ]
        
        tableGoGenesSplit <- split(tableGoGenes, tableGoGenes$ONTOLOGY) # Split by ontology
        tableGoGenesSplit <- tableGoGenesSplit[which(paste("GO_", names(tableGoGenesSplit), sep="") %in% annotations)] # Select the requested ontologies
        
        geneID2GO <- lapply(tableGoGenesSplit, function(x) split(as.character(x$GO),x[,geneIDtype]))
        names(geneID2GO) <- paste("GO_", names(geneID2GO), sep="")
        GO2geneID <- lapply(tableGoGenesSplit, function(x) split(x[,geneIDtype], as.character(x$GO)))
        names(GO2geneID) <- paste("GO_", names(GO2geneID), sep="")

        ret$genes2terms <- c(ret$genes2terms, geneID2GO)
        ret$terms2genes <- c(ret$terms2genes, GO2geneID)
        ret$asTable <- c(ret$asTable, list(GO=tableGoGenes))
     
        if(loadInstPkg("GO.db")) 
        {
            termsLabels <- AnnotationDbi::select(GO.db, keys=unlist(sapply(ret$terms2genes, names)), columns="TERM", keytype="GOID")  # ONLY annotations NEEDS TO BE GO!
            termsLabels <- split(termsLabels$TERM, termsLabels$GOID)# (as list for consistency with other annotationss)
            termsLabels <- lapply(termsLabels, capitalize)
            
            if(termAsLabel)
            {
                # Concatenate the term ID with the term label
                ret$terms2genes <- lapply(ret$terms2genes, function(x){
                    names(x) <- paste(names(x), ":", termsLabels[names(x)], sep="")
                    x})                        
            }else{ 
                # Add a field with the term name
                if(!"termID2name" %in% names(ret)) ret <- c(ret, list(termID2name=NULL))
                ret$termID2name <- c(ret$termID2name, list(GO=termsLabels))
            }
        }else{
            if(termAsLabel) warning("To add terms as gene set label GO.db package is required.")
        }
    }
    
    #### Build KEGG tables:
    if("KEGG" %in% annotations)
    {
        if("PATH" %in% AnnotationDbi::columns(org.db))
        {
            tableKeggGenes <- suppressWarnings(AnnotationDbi::select(org.db, keys=allGenes, columns="PATH", keytype=geneIDtype))
            tableKeggGenes <- tableKeggGenes[!is.na(tableKeggGenes[,"PATH"]),] # Bueno o malo?
            
            genes2kegg <- split(as.character(tableKeggGenes$PATH),tableKeggGenes[,geneIDtype])
            kegg2genes <- split(tableKeggGenes[,geneIDtype], as.character(tableKeggGenes$PATH))
    
            ret$genes2terms <- c(ret$genes2terms, list(KEGG=genes2kegg))
            ret$terms2genes <- c(ret$terms2genes, list(KEGG=kegg2genes))
            ret$asTable <- c(ret$asTable, list(KEGG=tableKeggGenes))
            
            if(loadInstPkg("KEGG.db")) 
            {
                # Labels...
                termsLabels <- AnnotationDbi::as.list(KEGG.db::KEGGPATHID2NAME)
                termsLabels <- termsLabels[which(names(termsLabels) %in% names(kegg2genes))]
                # termsLabels <- sapply(termsLabels, function(x) x[1]) # Take only first
                
                if(termAsLabel)
                {
                    # Concatenate the term ID with the term label
                    names(ret$terms2genes$KEGG) <- paste(names(ret$terms2genes$KEGG), ":", termsLabels[names(ret$terms2genes$KEGG)], sep="")        
                }else{ 
                    # Add a field with the term name
                    if(!"termID2name" %in% names(ret)) ret <- c(ret, list(termID2name=NULL))
                    ret$termID2name <- c(ret$termID2name, list(KEGG=termsLabels))
                } 
            }else{
                if(termAsLabel) warning("To add terms as gene set label KEGG.db package is required.")
            }
        }else{
            warning(paste("Kegg (pathway) info is not available in ", organismDb, sep=""))
        }
    }
    
    
    if("REACTOME" %in% annotations)
    {
        reactPw2geneID <- AnnotationDbi::as.list(reactome.db::reactomePATHID2EXTID)  
        geneID2reactPw <- AnnotationDbi::as.list(reactome.db::reactomeEXTID2PATHID)
        entrezIdsOrg <- AnnotationDbi::keys(org.db, keytype="ENTREZID")
        
        geneID2reactPw <- geneID2reactPw[which(names(geneID2reactPw) %in% entrezIdsOrg)]
        reactPw2geneID <- sapply(reactPw2geneID, function(x) x[which(x %in% entrezIdsOrg)])
        reactPw2geneID <-reactPw2geneID[sapply(reactPw2geneID, length)>0]
        
        tableReactGenes <- unlist(mapply(rep, names(reactPw2geneID), sapply(reactPw2geneID, length))) # equivalent to melt(reactPw2geneID)
        tableReactGenes <- data.frame(cbind(unlist(reactPw2geneID), tableReactGenes))
        colnames(tableReactGenes) <- c("ENTREZID", "PATH")
        
        
        # TRANSLATE?
        # TO DO!!
        if(geneIDtype != "ENTREZID") warning("At the moment REACTOME gene sets can only be provided for ENTREZID.")
        
        ret$genes2terms <- c(ret$genes2terms, list(REACTOME=geneID2reactPw))
        ret$terms2genes <- c(ret$terms2genes, list(REACTOME=reactPw2geneID))
        ret$asTable <- c(ret$asTable, list(REACTOME=tableReactGenes))
        
        # Labels...
        termsLabels <- AnnotationDbi::as.list(reactomePATHID2NAME)
        termsLabels <- termsLabels[which(names(termsLabels) %in% names(reactPw2geneID))]
        termsLabels <- sapply(termsLabels, function(x) x[1]) # Take only first
        termsLabels <- strsplit(termsLabels, ": ")
        termsLabels <- sapply(termsLabels, function(x){
            paste(paste(x[-1],collapse=" ")," (",x[1],")",sep="")
        })
        if(termAsLabel)
        {
            # Concatenate the term ID with the term label
            names(ret$terms2genes$REACTOME) <- paste(names(ret$terms2genes$REACTOME), ":", termsLabels[names(ret$terms2genes$REACTOME)], sep="")        
        }else{ 
            # Add a field with the term name
            if(!"termID2name" %in% names(ret)) ret <- c(ret, list(termID2name=NULL))
            ret$termID2name <- c(ret$termID2name, list(REACTOME=termsLabels))
        }        
    }    
    return(ret)
}



# # termList <- sigTerm$goid
# # organism <- "org.Hs.eg.db"
# # geneIDtype<-"ENTREZID"
# # geneList <- entrezIDch
# 
# ### Valido para otros FEA/GSEA cuando no dan genes en el gtset... etc
# findGenesInTerm <- function(termIds, geneIds, geneIDtype, organismDb, addTermLabel=TRUE)
# {
#     termIds <- as.character(termIds)
#     dbs <- buildDatabases(organismDb,geneIDtype)
#     
#     ret <- cbind(TermID=termIds, Genes=sapply(termIds, function(term) 
#     {
#         genes <- unique(c(dbs$tableGoGenes[which(dbs$tableGoGenes$GO==term),geneIDtype], # GO
#                    dbs$tableKeggGenes[which(dbs$tableKeggGenes$PATH==term),geneIDtype])) # KEGG
#         return(paste(as.character(genes[genes %in% geneIds]), collapse=","))       
#     }))
#     
#     
#     if(addTermLabel && ("goTerms" %in% names(dbs)))
#     {
#         termLabel <- dbs$goTerms[ret[,"TermID"],]
#         termLabel <- apply(termLabel, 1, function(x) paste(x[1],"~", x[2], sep=""))
#         ret <- cbind(TermID=ret[,"TermID"], Term=termLabel, Genes=ret[,"Genes"])
#     }else{
#         colnames(ret) <- c("Term", "Genes")
#     }
#         
#     return(ret)
# }  


