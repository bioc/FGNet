# also used in toMatrix & clDescriptions
formatTerms<- function(termTbl)
{
    if(any(grepl(":",termTbl, fixed=TRUE)))
    {
        termTbl <- strsplit(termTbl, ":")
        
        ret <- lapply(termTbl, function(x) 
        {
            if(length(x)<3)
            {
                capitalize(paste(paste(x[length(x)],sep=""))) # Drops annotation
            }else# format ANNOT:ID:Term description...
            {
                capitalize(paste(paste(x[3:length(x)], collapse=" "), " (",paste(x[1:2], collapse=":"),")",sep=""))
            }
        })
    }else
    {
        ret <- termTbl
    }
    return(ret)
}

# keepTermsCol: Column name. Only geneTerm sets with this column TRUE will be kept.  (Used by gtLinker and GAGE. Can also be used to filter...)
clustersTable <- function(geneTermSets, clusterColumn="Cluster", colsAvg=NULL, keepTermsCol=NULL, addKeyWordsTerm=TRUE, sortBy=NULL, decreasing=TRUE)
{
    ###########################
    # Create clusters table
    tablaClusters <- NULL 
    
    if(clusterColumn %in% names(geneTermSets))
    {
        if(!is.null(keepTermsCol)) geneTermSets <- geneTermSets[geneTermSets[,keepTermsCol],,drop=FALSE]
        
        for(cl in unique(geneTermSets[,clusterColumn]))
        {
            tmpTable <- as.matrix(geneTermSets[which(geneTermSets[,clusterColumn] == cl),c(clusterColumn, "Terms", "Genes", colsAvg)])
    
            # Collapse Terms
            tmpTerms <- tmpTable[,"Terms", drop=FALSE]
            if(addKeyWordsTerm) keyWordsTerm <- capitalize(keywordsTerm(list(cbind(TermDescription=sapply(strsplit(tmpTerms, ":"), function(x) x[length(x)]))), nChar=100))
            tmpTerms <- formatTerms(tmpTerms)
            tmpTerms <- paste(tmpTerms, collapse=";")
    
            # Merge genes from geneTerm sets
            tmpGenes <- list()
            for(i in 1:length("Genes"))
            {
                tmpGenes[[i]] <- sort(unique(unlist(strsplit(as.character(tmpTable[,"Genes"[i]]), split=","))))
                nGenes <- length(tmpGenes[[i]])
                tmpGenes[[i]] <- paste(tmpGenes[[i]], collapse=",")
            }
            tmpGenes <- unlist(tmpGenes)
            
            # Calculate average from requested colums
            colsAvgCluster <- NULL
            if(!is.null(colsAvg)) 
            {
                colsAvgCluster <- sapply(colsAvg, function(x) 
                    {
                        x <- tmpTable[,x]
                        if(all(!is.na(suppressWarnings(as.numeric(x))))) 
                        {
                            x <- mean(as.numeric(x))
                        }else
                        {
                            paste(unique(x), collapse=", ")
                        }
                        
                    })
            }
            rowCluster <- c(Cluster=cl, nGenes=nGenes, colsAvgCluster, Genes=tmpGenes, Terms=tmpTerms)
            if(addKeyWordsTerm) rowCluster <- c(rowCluster, keyWordsTerm=keyWordsTerm)
            tablaClusters <- rbind(tablaClusters, rowCluster)
        }
        colnames(tablaClusters)[1] <- clusterColumn
        rownames(tablaClusters) <- as.character(tablaClusters[,clusterColumn])
        tablaClusters <- data.frame(tablaClusters)
        
        
        if(!is.null(sortBy))
        {
            sortBy <- as.character(tablaClusters[,sortBy])
            if(all(!is.na(suppressWarnings(as.numeric(sortBy))))) sortBy <- as.numeric(sortBy)
            tablaClusters <- tablaClusters[order(sortBy,decreasing=decreasing),]
        }
    }
    return(tablaClusters)        
}