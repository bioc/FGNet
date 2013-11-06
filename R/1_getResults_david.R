
#results <- getResults_David (fileName)
getResults_david <- function(inputFileLocation, path=getwd(), jobName="")
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
			file.copy(inputFileLocation, downloadedFileName)
			message(paste(inputFileLocation, " copied to ", downloadedFileName, sep=""))
			inputFileLocation <- downloadedFileName
		}
	}
	
	# Read file & process
	if (!file.exists(inputFileLocation))
	{
		stop("Can't open or download the results file.")
	}
	else
	{
		inputFile <- file(inputFileLocation, "rt")
		lineas <- readLines(inputFile) 
	
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
			linea
		}

		close(inputFile) 
		tablaGeneTermSets <- data.frame(tablaGeneTermSets)
		
	
		# Create clusters table
		
		tablaClusters <- NULL
		for(cluster in names(clusterScore))
		{
			tmpTable <- as.matrix(tablaGeneTermSets[which(tablaGeneTermSets[,"Cluster"] == cluster),c("Cluster", "Term", "Category", "Genes")])
			# Kegg
			keggs <- which(tmpTable[,"Category"] =="KEGG_PATHWAY")
			tmpTable[keggs, "Term"] <- paste("KEGG:", tmpTable[keggs, "Term"], sep="")				# Used in createHtml
				
			# Not Kegg or GO
			gos <- grep("GOTERM", tmpTable[,"Category"])
			iprs <- grep("INTERPRO", tmpTable[,"Category"])
			otherAnnot <- which(!(1:dim(tmpTable)[1] %in% c(keggs, gos, iprs)))
			tmpTable[otherAnnot, "Term"] <- sub(":", ". ", tmpTable[otherAnnot, "Term"])
			tmpTable[otherAnnot, "Term"] <- sub(";", ". ", tmpTable[otherAnnot, "Term"])
			tmpTable[otherAnnot, "Term"] <- paste(tmpTable[otherAnnot, "Category"], tmpTable[otherAnnot, "Term"], sep=":")	
						
			# All
			tmpTerms <- paste(sub("~", ":", tmpTable[,"Term"]), collapse=";")
			tmpGenes <- unique(unlist(strsplit(as.character(tmpTable[,"Genes"]), split=", ")))
			nGenes <- length(tmpGenes)
			tmpGenes <- paste(tmpGenes, collapse=",")
			tablaClusters <- rbind(tablaClusters, c(cluster, clusterScore[cluster], nGenes, tmpGenes, tmpTerms))
		}
		colnames(tablaClusters) <- c("Cluster", "EnrichmentScore", "nGenes", "Genes", "Terms")
		tablaClusters <- data.frame(Cluster=tablaClusters[,"Cluster"], EnrichmentScore=as.numeric(tablaClusters[,"EnrichmentScore"]), nGenes=as.numeric(tablaClusters[,"nGenes"]), Genes=tablaClusters[,"Genes"], Terms=tablaClusters[,"Terms"], stringsAsFactors = FALSE)
	}
	return(list(clusters=tablaClusters, geneTermSets=tablaGeneTermSets, fileName=inputFileLocation))
}
