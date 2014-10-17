#filterAttribute: Sorted by metagroup/Cluster!!

fea2incidMat <- function(feaResults, key="Genes", sepChar=NULL, clusterColumn=NULL, filterAttribute=NULL, filterOperator="<", filterThreshold=0, removeFilteredGtl=NULL)
{
    # Check arguments
    grType <- "Cl"
    metagroups <- NULL
    if(is.list(feaResults) && ("metagroups" %in% names(feaResults))) grType <- "Mg"
    if(is.list(feaResults) && ("metagroups" %in% names(feaResults))) metagroups <- feaResults$metagroups   # Para david no hace falta, son lo mismo
    if(is.list(feaResults) && ("geneTermSets" %in% names(feaResults))) geneTermSets <- feaResults$geneTermSets

    if(is.data.frame(feaResults) || is.matrix(feaResults)) 
    {
        geneTermSets <- data.frame(feaResults)
    }
    
    # Fom other types of tools: Consider each gtset an individual cluster (assign each of them a color)
    if(!any(c("Cluster", "Metagroup") %in% colnames(geneTermSets)))
    {
        geneTermSets <- data.frame(cbind(Gtset=1:nrow(geneTermSets),geneTermSets))
        clusterColumn <- "Gtset"
    }
    
    if(!is.data.frame(geneTermSets)) stop("feaResults should be the formatted gene-term sets.")
    if(dim(geneTermSets)[1] == 0) stop("0 gene term sets.")
    if(!is.na(suppressWarnings(as.numeric(filterThreshold)))) filterThreshold <- as.numeric(filterThreshold)
    #if(!is.numeric(filterThreshold)) stop("filterThreshold should be a number.")
    if(!is.null(filterAttribute) && is.na(filterAttribute)) 
    {
        filterAttribute <- NULL
        filterOperator <- NULL 
        filterThreshold <- 0        
    }
    if(is.null(filterThreshold)) filterThreshold <- 0 
    if(!is.null(filterAttribute) && !is.data.frame(filterAttribute)) 
    {
        if(class(feaResults)=="list")
        {
            tmpClusterResults <- feaResults[[which(names(feaResults) %in% c("clusters","metagroups"))]]
            if(is.null(tmpClusterResults)) tmpClusterResults <- geneTermSets
        }else
        {
            tmpClusterResults <- geneTermSets
        }        
        if(is.character(filterAttribute) && (filterAttribute %in% colnames(tmpClusterResults)))
        {
            tmpFilterAttribute <- tmpClusterResults[,filterAttribute, drop=FALSE]  
            tmpFilterAttribute <- setNames(as.character(tmpFilterAttribute[,1]), rownames(tmpFilterAttribute))
            if(!any(is.na(suppressWarnings(as.numeric(tmpFilterAttribute))))) tmpFilterAttribute <- setNames(as.numeric(tmpFilterAttribute), names(tmpFilterAttribute))
            
            tmpFilterAttribute <- data.frame(tmpFilterAttribute)
            colnames(tmpFilterAttribute) <- filterAttribute
            filterAttribute <- tmpFilterAttribute
            
        } else stop("filterAttribute should be a column from the data.frame returned by one of the fea_ functions.")    
    }
    if(is.null(filterAttribute) && filterThreshold!=0) stop("To filter provide an filterAttribute.")
    key <- capitalize(tolower(key))
    
    if(is.null(clusterColumn))
    {
        clusterColumn <- "Cluster"
        if(("metagroups" %in% names(feaResults)) && ("Metagroup" %in% colnames(feaResults$geneTermSets))) clusterColumn <- "Metagroup"
    }
    
    #if(!key %in% c("Genes", "Terms")) stop("key should be either 'Genes' or 'Terms'.")
    if(is.null(removeFilteredGtl)) removeFilteredGtl <- TRUE
    if(!is.logical(removeFilteredGtl)) stop("removeFilteredGtl should be either TRUE or FALSE.")
    if(removeFilteredGtl && is.null(metagroups))  removeFilteredGtl <- FALSE
    if(removeFilteredGtl && key!="Terms") removeFilteredGtl <- FALSE
      
    # Initialize
    if(is.null(sepChar))
    {
        if(key == "Genes") sepChar <- ","
        if(key == "Terms") sepChar <- ";"
    }
    
    emptyGtSet<-nchar(as.character(geneTermSets[,key]))==0 # Any gtset does not have genes / terms?
    if(any(emptyGtSet))
    {
        warning("There are gene-term sets withouth genes.")
        geneTermSets <- geneTermSets[!emptyGtSet,]
    }    

    colsNeeded <- which(colnames(geneTermSets) %in% c("Cluster", "Metagroup", key, clusterColumn))
    geneTermSets <- as.matrix(geneTermSets[,colsNeeded, drop=FALSE]) # Group & Genes
    if(all(!is.na(suppressWarnings(as.numeric(geneTermSets[,clusterColumn]))))) geneTermSets[,clusterColumn] <- gsub(" ", "", geneTermSets[,clusterColumn])

    # Filter & sort by filterAttribute    
    filtrar <- NULL
    if(!is.null(filterAttribute))
    {  
        # Filter
        if(is.character(filterThreshold)) filterThreshold <- paste("\"",filterThreshold,"\"",sep="")

        filtrar <- rownames(filterAttribute)[eval(parse(text=paste("which(filterAttribute[,1]", filterOperator, filterThreshold,")",sep=""))) ]
        if(length(filtrar)>0)
        {
            # Filter by Cluster or Gene-Term set?
            filtByCluster <- sum(filtrar %in% geneTermSets[,clusterColumn]) >= sum(filtrar %in% rownames(geneTermSets))
            if(filtByCluster)  
            {
                 geneTermSets <- geneTermSets[-which(geneTermSets[,clusterColumn] %in% filtrar),, drop=FALSE]
                  message(paste("The following metagroups/clusters have been filtered out: ", paste(filtrar, collapse=", "), sep=""))
            }else{     
                 geneTermSets <- geneTermSets[-which(rownames(geneTermSets) %in% filtrar),, drop=FALSE]
                 message(paste("The following gene-term sets have been filtered out: ", paste(filtrar, collapse=", "), sep=""))
            }
            filterAttribute <- filterAttribute[-which(rownames(filterAttribute) %in% filtrar), , drop=FALSE]
        
            if(nrow(geneTermSets)==0) stop(paste("All metagroups/clusters have ", colnames(filterAttribute), filterOperator, filterThreshold,". Try modifying the filter threshold.", sep=""))
            
            # Sort
            # filterAttribute <- filterAttribute[order(filterAttribute[,1], decreasing=TRUE), , drop=FALSE]
            
            if(filtByCluster)
            {
                sortedGeneTermSets <- NULL
                for(gr in rownames(filterAttribute)) sortedGeneTermSets <- rbind(sortedGeneTermSets, geneTermSets[ geneTermSets[,clusterColumn] %in% gr,,drop=FALSE])
                geneTermSets <- sortedGeneTermSets
            }else
            {
                geneTermSets <- geneTermSets[rownames(filterAttribute),]
            }
            # Filter MG:
            if(!is.null(metagroups)) metagroups <- metagroups[which(!rownames(metagroups) %in% filtrar),]
        }
    } else
    {
      ord <- geneTermSets[,clusterColumn]
      if(all(!is.na(suppressWarnings(as.numeric(ord))))) ord <-as.numeric(ord)
      geneTermSets <- geneTermSets[order(ord),, drop=FALSE]        
    }

    # Clear (for DAVID)
    # Remove spaces in genes
    if(key=="Genes") geneTermSets[,key] <- sapply(geneTermSets[,key], function(tmp) gsub(pattern=" ", replacement="",  x=tmp, fixed=TRUE))
    # Change ~ GO separator:
    if(key=="Terms") geneTermSets[,key] <- sapply(geneTermSets[,key], function(tmp) gsub(pattern="~", replacement=":",  x=tmp, fixed=TRUE))
    
    ################################
    # Transform into tables
    ################################
                                
    ### GeneTermSets
    # List of genes/terms in each Gene-Term Set
    gtSetsList <- as.list(as.character(geneTermSets[,key]))
    if(!is.null(sepChar)) 
    {
        gtSetsList <- lapply(gtSetsList, function(gr) { 
                unlist(strsplit(as.character(gr), sepChar, fixed=TRUE))
        })
       
    }    
                         
    # Columna 1: Metagrupo/Cluster al que pertenece
    names(gtSetsList) <- paste(as.character(geneTermSets[,clusterColumn]), ".", unlist(lapply(table(factor(geneTermSets[,clusterColumn]))[unique(geneTermSets[,clusterColumn])], function(x) 1:x)), sep="")
     
    # Remove filtered terms from metagroups (only GeneTermLinker)
    if(removeFilteredGtl && key=="Terms")
    {
        nonFilteredList <- getMGTerms(metagroups, grType="Metagroup")
        nonFilteredList <- sapply(nonFilteredList, function(gr){
          paste(gr[,"TermID"], ":", gr[,"TermDescription"], sep="")
        })
        names(nonFilteredList) <- gsub("Metagroup ", "", names(nonFilteredList))
        
        newGtSetsList <- NULL
        for(mg in names(nonFilteredList))
        {
            gtsets <- grep(paste("^",mg, "[.]", sep=""), names(gtSetsList))
            gtSetsList[gtsets] <- lapply(gtSetsList[gtsets], function(gtset) gtset[tolower(gtset) %in% tolower(nonFilteredList[[mg]])])
            tmp <- unique(gtSetsList[gtsets])
            names(tmp) <- paste(mg, ".", 1:length(tmp), sep="")
            newGtSetsList <- c(newGtSetsList, tmp)
        }
        gtSetsList <- newGtSetsList
    }
    if(key=="Terms") gtSetsList <- sapply(gtSetsList, function(x) unlist(formatTerms(x)))
    
    nKeys <- length(unique(unlist(gtSetsList)))

    # Genes - GeneTermSet matrix
    nGtSets <- length(gtSetsList)
    gtSetsMatrix <- matrix(ncol=nGtSets, nrow=nKeys, data=0)
    colnames(gtSetsMatrix) <- names(gtSetsList)
    rownames(gtSetsMatrix) <- sort(unique(unlist(gtSetsList)))
    for( gts in  names(gtSetsList) )
    {
        gtSetsMatrix[eval(gtSetsList[[gts]]),gts]<-1
    }

    # Metagroups / Clusters     (Based on GeneTermSets. In gtLinker it might include filtered terms)
    mgList <- list()
    for(mg in unique(geneTermSets[,clusterColumn])) # Non-filtered mg
    {
        mgList[[mg]] <- unlist(gtSetsList[which(sapply(strsplit(names(gtSetsList),".",fixed=TRUE), function(x) x[-length(x)]) == mg)])
        if(!is.null(sepChar))
        { 
            mgList[[mg]] <- unique(unlist(lapply(mgList[[mg]], function(mg) { strsplit(as.character(mg), sepChar, fixed=TRUE)}))) 
        }else
        {
            mgList[[mg]] <- unique(unlist(mgList[[mg]])) 
        }

        if(key=="Genes") mgList[[mg]] <- gsub(pattern=" ", replacement="",  x=mgList[[mg]], fixed=TRUE)
    }
    
    nGroups <- length(mgList)
    nKeys <- length(unique(unlist(mgList)))
    metagroupsMatrix <- matrix(ncol=nGroups, nrow=nKeys, data=0)
    colnames(metagroupsMatrix) <- names(mgList)
    rownames(metagroupsMatrix) <- sort(unique(unlist(mgList)))
    for( mg in  names(mgList) )
    {
        metagroupsMatrix[eval(mgList[[mg]]),mg]<-1
    }

    ret <- list(metagroupsMatrix=metagroupsMatrix, gtSetsMatrix=gtSetsMatrix, filteredOut=filtrar)
    if(grType=="Cl") names(ret)[1] <- "clustersMatrix"
    return(ret)
}

