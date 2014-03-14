
#results <- getResults_David (fileName)
getResults_david <- function(inputFileLocation, path=getwd(), jobName="", geneLabels=NULL)
{
	# Check arguments
	if(!is.character(inputFileLocation)) stop("inputFileLocation not valid.")
	if(!file.exists(path)) stop("The given path does not exist.")
	if(!is.character(jobName)) stop("jobName should be character.")
	
	# Create folder
	if(substring(path, first=nchar(path), last=nchar(path)) != "/") path <- paste(path, "/", sep="")
	folder <- path
	if(jobName != "") 
	{
		folder <- paste(path, jobName, "/", sep="")	
		jobName <- paste(jobName, "_", sep="")
	}

	if ((!file.exists(folder)) && (!folder=="")){
		dir.create(file.path(folder))		
	}
	
	# Download file
	if(grepl("http://", inputFileLocation)) 
	{
		downloadedFileName <- paste(folder, jobName,"DavidClustering.txt", sep="")
		if(!url.exists(inputFileLocation)) stop("Download URL (inputFileLocation) is not available.")
		download.file(inputFileLocation, destfile=downloadedFileName, quiet=TRUE)
		inputFileLocation <- downloadedFileName
		message(paste("Raw results downloaded and saved as ", inputFileLocation, sep=""))
	}else{
		if(jobName!= "")
		{
			downloadedFileName <- paste(folder, jobName, "DavidClustering.txt", sep="")
			file.copy(inputFileLocation, downloadedFileName, overwrite=TRUE)
			message(paste(inputFileLocation, " copied to ", downloadedFileName, sep=""))
			inputFileLocation <- downloadedFileName
		}
	}

	# Read file & process
	if (!file.exists(inputFileLocation))
	{
		stop("Can't open or download the results file.")
	}	else
	{
		inputFile <- file(inputFileLocation, "rt")
			lineas <- readLines(inputFile)
			if(length(lineas)==0) stop("DAVID returned 0 clusters.")
		
			columns <- c("Cluster", strsplit(lineas[2], "\t", fixed=TRUE)[[1]])
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
					tablaGeneTermSets <- rbind(tablaGeneTermSets, c(paste(cluster, sep=""), linea))
				}
			}
		close(inputFile) 
		
		# Replace Gene Names?
		if(!is.null(geneLabels))
		{
			if(length(geneLabels) != length(unique(geneLabels))) stop("geneLabels IDs are not unique.")
			if(length(geneLabels) != length(unique(names(geneLabels)))) stop("geneLabels names are not unique.")
			
			colGenes <- which(colnames(tablaGeneTermSets) == "Genes")
			tablaGeneTermSets <- tablaGeneTermSets[,c(colnames(tablaGeneTermSets)[-colGenes], "Genes")]
			tablaGeneTermSets <- cbind(tablaGeneTermSets, GenesIDs=tablaGeneTermSets[,"Genes"])
			for( i in 1:length(geneLabels))
			{
				tablaGeneTermSets[,"Genes"] <- sapply(tablaGeneTermSets[,"Genes"], function(x) sub(geneLabels[i], names(geneLabels[i]), x))
			}
		}
		tablaGeneTermSets <- data.frame(tablaGeneTermSets)
		
		# Create clusters table
		tablaClusters <- NULL
		# Gene labels & IDS?
		colGenes <- "Genes"
		if(!is.null(geneLabels)) colGenes <- c(colGenes, "GenesIDs")
		
		for(cluster in names(clusterScore))
		{
			tmpTable <- as.matrix(tablaGeneTermSets[which(tablaGeneTermSets[,"Cluster"] == cluster),c("Cluster", "Term", "Category", colGenes)])
			# Kegg
			keggs <- which(tmpTable[,"Category"] == "KEGG_PATHWAY")
			tmpTable[keggs, "Term"] <- paste("KEGG:", tmpTable[keggs, "Term"], sep="")				# Used in createHtml
			
			# Reactome
			reactomes <- grep("REACTOME_PATHWAY", tmpTable[,"Category"]) 
			tmpTable[reactomes, "Term"] <- sub("REACT_", "REACT:", tmpTable[reactomes, "Term"])
			# Not Kegg or GO or...
			gos <- grep("GOTERM", tmpTable[,"Category"])
			iprs <- grep("INTERPRO", tmpTable[,"Category"])		
			
			otherAnnot <- which(!(1:dim(tmpTable)[1] %in% c(keggs, gos, iprs, reactomes)))
			if (length(otherAnnot) > 0)
			{
				#tmpTable[otherAnnot, "Term"] <- gsub(":", ". ", tmpTable[otherAnnot, "Term"]) -> Pasar el codigo al final con parentesis
				tmpTable[otherAnnot, "Term"] <- sapply(strsplit(tmpTable[otherAnnot, "Term"], split=":"), function(x) paste(x[-1], " (ID ", x[1], ")", sep=""))
				tmpTable[otherAnnot, "Term"] <- paste(tmpTable[otherAnnot, "Category"], tmpTable[otherAnnot, "Term"], sep=":")	
			}
			tmpTable[, "Term"] <- gsub(";", ". ", tmpTable[, "Term"])
			
			# Genes & terms
			tmpTerms <- paste(sub("~", ":", tmpTable[,"Term"]), collapse=";")
			tmpGenes <- list()
			for(i in 1:length(colGenes))
			{
				tmpGenes[[i]] <- unique(unlist(strsplit(as.character(tmpTable[,colGenes[i]]), split=", ")))
				nGenes <- length(tmpGenes[[i]])
				tmpGenes[[i]] <- paste(tmpGenes[[i]], collapse=",")
			}
			tmpGenes <- unlist(tmpGenes)
			
			tablaClusters <- rbind(tablaClusters, c(cluster, clusterScore[cluster], nGenes, tmpGenes, tmpTerms))
		}
		colnames(tablaClusters) <- c("Cluster", "EnrichmentScore", "nGenes", colGenes, "Terms")
		tablaClusters <- data.frame(Cluster=tablaClusters[,"Cluster"], EnrichmentScore=as.numeric(tablaClusters[,"EnrichmentScore"]), nGenes=as.numeric(tablaClusters[,"nGenes"]), tablaClusters[,c(colGenes, "Terms"), drop=FALSE], stringsAsFactors = FALSE)
		rownames(tablaClusters) <- tablaClusters[,"Cluster"]
	}
	
	return(list(clusters=tablaClusters, geneTermSets=tablaGeneTermSets, fileName=inputFileLocation))
}
