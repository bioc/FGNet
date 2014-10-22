# results <- getResults_David (fileName)
# copyFile (for internal use)
format_david <- function(fileName, jobName=NULL, geneLabels=NULL, moveFile=FALSE)
{
    # Check arguments
    if(!is.character(fileName)) stop("fileName not valid.")
    if(!is.null(jobName) && !is.character(jobName)) stop("jobName should be character.")
    
    ###############################
    # Prepare jobName / folder
    isURL <- grepl("http://", fileName)
    
    if(is.null(jobName)) 
    {
        if(isURL) jobName <- paste(sample(100000:999999,size=1), "_DavidClustering",sep="")
        if(!isURL) jobName <- sub("_raw", "", getJobName(fileName), fixed=TRUE)
    }
    folder <- jobName

    # Create folder
    if((!file.exists(file.path(folder))))
    {
        dir.create(file.path(folder))        
    }
    currWD <- getwd()
    setwd(folder)
    
    ###############################
    # Download file    
    if(isURL) 
    {
        setwd(currWD)
        downloadedFileName <- paste(jobName, "_raw.txt", sep="")
        download.file(fileName, destfile=downloadedFileName, quiet=TRUE)
        fileName <- downloadedFileName
        message(paste("Raw results downloaded and saved as ", fileName,"'", sep=""))
    }else{
        if(!file.exists(fileName)) 
        {
            setwd(currWD)
            stop("The file does not exist.")
        }
        downloadedFileName <- paste(jobName, "_raw.txt", sep="")
        if(fileName != paste(getwd(),downloadedFileName, sep=.Platform$file.sep)) # Same file, no need to copy
        {
            tmp <- paste(sample(1:9999999,1), ".txt", sep="") # Needed in case origin and destination are the same file (i.e. link)
            file.copy(fileName, tmp, overwrite=TRUE) 
            file.copy(tmp, downloadedFileName, overwrite=TRUE) 
            file.remove(tmp)
            if(moveFile) file.remove(fileName) 
            message(paste(fileName, " copied to '", folder,.Platform$file.sep, downloadedFileName,"'", sep=""))
        }
        fileName <- downloadedFileName
    }
    
    # Read file & process
    if (!file.exists(fileName))
    {
        setwd(currWD)
        stop("Can't open or download the results file.")
    }    else
    {
        inputFile <- file(fileName, "rt")
        lineas <- readLines(inputFile)
        if(length(lineas)==0) stop("DAVID returned 0 clusters.")
        
        columns <- c("Cluster", "ClusterEnrichmentScore", strsplit(lineas[2], "\t", fixed=TRUE)[[1]])
        tablaGeneTermSets <- matrix(NA, ncol=length(columns), nrow=0, dimnames=list(c(), columns))
        clusterScore <- NULL
        cluster <- 0
        
        lineas <- lineas[which(lineas != "")]
        lineas <- lineas[which(lineas != lineas[2])]
        lineas <- strsplit(lineas, "\t", fixed=TRUE)
        
        for(linea in lineas)
        {
            if(length(linea)==2)
            {
                cluster <- cluster+1
                clusterScore[[cluster]] <- as.numeric(strsplit(linea[2], ": ", fixed="TRUE")[[1]][2])
                names(clusterScore)[cluster] <- cluster
            }
            else 
            {
                tablaGeneTermSets <- rbind(tablaGeneTermSets, c(cluster, clusterScore[cluster], linea))
            }
        }
        close(inputFile) 
        tablaGeneTermSets <- data.frame(tablaGeneTermSets)
        
        ##########################
        # Format tablaGeneTermSets
        colnames(tablaGeneTermSets)[which(colnames(tablaGeneTermSets)=="Term")] <- "Terms" 
        tablaGeneTermSets$Terms <- sapply(tablaGeneTermSets$Terms, function(x) gsub("~",":",x))
        tablaGeneTermSets$Genes <- sapply(tablaGeneTermSets$Genes, function(x) gsub(" ","",x))
        
        # Kegg
        keggs <- which(tablaGeneTermSets[,"Category"] == "KEGG_PATHWAY")
        if(length(keggs)>0) tablaGeneTermSets[keggs, "Terms"] <- paste("KEGG:", tablaGeneTermSets[keggs, "Terms"], sep="")        
        
        # Reactome
        reactomes <- grep("REACTOME_PATHWAY", tablaGeneTermSets[,"Category"]) 
        if(length(reactomes)>0)  tablaGeneTermSets[reactomes, "Terms"] <- sub("REACT_", "REACT:", tablaGeneTermSets[reactomes, "Terms"])
        
        # Interpro
        iprs <- grep("INTERPRO", tablaGeneTermSets[,"Category"])
        if(length(iprs)>0) tablaGeneTermSets[iprs, "Terms"] <- sub("IPR", "IPR:", tablaGeneTermSets[iprs, "Terms"])
        
        # Not any of these annotations...
        gos <- grep("GOTERM", tablaGeneTermSets[,"Category"])        
        otherAnnot <- which(!(1:dim(tablaGeneTermSets)[1] %in% c(keggs, gos, iprs, reactomes)))
        if (length(otherAnnot) > 0)
        {
            tablaGeneTermSets[otherAnnot, "Terms"] <- sapply(strsplit(tablaGeneTermSets[otherAnnot, "Terms"], split=":"), function(x) paste(x[-1], " (ID ", x[1], ")", sep=""))
            tablaGeneTermSets[otherAnnot, "Terms"] <- paste(tablaGeneTermSets[otherAnnot, "Category"], tablaGeneTermSets[otherAnnot, "Terms"], sep=":")    
        }
        tablaGeneTermSets[, "Terms"] <- gsub(";", ". ", tablaGeneTermSets[, "Terms"])
        
        # Replace Gene Names?
        tablaGeneTermSets <- addGeneLabels(tablaGeneTermSets, geneLabels)
        
        # Save formatted file
        fileName <- paste(jobName, "_formatted.txt", sep="")
        write.table(tablaGeneTermSets, file=fileName, quote=FALSE, row.names=FALSE, sep="\t")  
    }
    ret <- readGeneTermSets(fileName, tool="DAVID")
    setwd(currWD)
    invisible(ret)
}
