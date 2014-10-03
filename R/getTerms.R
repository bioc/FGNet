# Return:
# $`Cluster 1`
# Terms                                      
# [1,] "Cellular component organization"          
# [2,] "Organelle organization"                   
# [3,] "Positive regulation of catalytic activity"
# [4,] "Positive regulation of molecular function"
# ...

getTerms <- function(feaResults, returnValue="description") #GO, KEGG
{    
    if(!is.list(feaResults)) stop("Please provide the whole output from getResults_david() or getResults_gtLinker().")  # Only Clusters/Metagroups are used. (no need of $GeneTermSets)
    
                            grNames <- NULL
                            if("clusters" %in% names(feaResults))
                            {    
                                if(!is.null(feaResults$clusters))    
                                {
                                    grType<-"Cluster"
                                    feaResults <- feaResults$clusters
                                    if("Cluster" %in% colnames(feaResults)) grNames <- feaResults$Cluster
                                }else{
                                    
                                    grType <- "Gene-term set"
                                    feaResults <- feaResults$geneTermSets
                                    feaResults$Terms <- unlist(formatTerms(as.character(feaResults$Terms)))
                                    
                                    grNames <- rownames(feaResults)
                                }
                            }
                            if("metagroups" %in% names(feaResults))
                            {
                                grType <- "Metagroup"
                                feaResults <- feaResults$metagroups
                                grNames <- rownames(feaResults)
                            }

    mgTerms <- getMGTerms(feaResults, grType=grType)
    
    if(returnValue=="GO") # Return go IDs
    {        
        terms <- sapply(mgTerms, function(x) x[grepl(pattern="GO:",x[,"TermID"]),"TermID"])
        #terms <- sapply(terms, function(x) sapply(x, function(y) strsplit(y,split=":[[:alpha:]]",perl=TRUE)[[1]][1]))
        names(terms) <- paste(grType, grNames)

    }else if(returnValue=="KEGG") # Return KEGG IDs
    {    
        terms <- sapply(mgTerms, function(x) 
            {
            x <- x[grepl(pattern="KEGG:",x[,"TermID"]),,drop=FALSE]
            x <- setNames(x[,"TermID"], x[,"TermDescription"])
            })
        terms <- sapply(terms, function(x) sapply(x, function(y) strsplit(y,split=":",perl=TRUE)[[1]][2]))
        names(terms) <- paste(grType, grNames)
        # Quitar vacios?
        
    }else if(returnValue=="description") # Return term description
    {    
        terms <- lapply(mgTerms, function(x) 
            {
            temp <- rbind(x)[,"TermDescription", drop=FALSE]
            colnames(temp) <- "Terms"
            rownames(temp) <- NULL
            return(temp)
            })
    }
    return(terms)    
}
