
#tables <- toTables_gtLinker( results$globalMetagroups, results$geneTermSets)
#tables <- adjMatrix(results$geneTermSets, results$globalMetagroups[,"Silhouette Width"])

#attribute: Sorted by metagroup/Cluster!!

adjMatrix <- function(geneTermSets, attribute=NULL, threshold=0)
{
	# Check arguments
	if(!is.data.frame(geneTermSets)) stop("geneTermSets should be the data.frame returned by getResults_gtLinker().")
	if(!is.numeric(threshold)) stop("threshold should be a number.")
	if(!is.null(attribute) && !is.data.frame(attribute)) stop("attribute should be the data.frame returned by getResults...().")	
		if(is.null(attribute) && threshold!=0) stop("To filter provide an attribute.")
	 	
	cols <- which(colnames(geneTermSets) %in% c("Cluster", "Metagroup", "Genes"))
	geneTermSets <- as.matrix(geneTermSets[,cols]) # Group & Genes
	
	# Filter & sort by attribute	
	filtrar <- NULL
	if(!is.null(attribute))
	{	
		# Filter
		filtrar <- rownames(attribute)[which(attribute[,1]<threshold)]
		
		if(length(filtrar)>0)
		{
		 	geneTermSets <- geneTermSets[-which(geneTermSets[,1] %in% filtrar),]
		 	attribute <- attribute[-which(rownames(attribute) %in% filtrar), , drop=FALSE]
		 	message(paste("The following metagroups/clusters have been filtered out: ", paste(filtrar, collapse=", "), sep=""))
		}
		if(dim(geneTermSets)[1]==0) stop(paste("There are no metagroups/clusters over the threshold (",threshold,"). Try lowering it.", sep=""))
		
		# Sort
		attribute <- attribute[order(attribute[,1], decreasing=TRUE), , drop=FALSE]
		sortedGeneTermSets <- NULL
		for(gr in rownames(attribute))
		{
			sortedGeneTermSets <- rbind(sortedGeneTermSets, geneTermSets[ geneTermSets[,1] %in% gr,])
		}
		geneTermSets <- sortedGeneTermSets
	}
	
	 
	################################
	# Transform into tables
	################################
								
	### GeneTermSets
	# Columna 1: Metagrupo/Cluster al que pertenece
	nGtSets<- dim(geneTermSets)[1]
	
	# List of genes in each Gene-Term Set
	gtsetGenesList <- as.character(geneTermSets[,"Genes"])
	#gtsetGenesList <- lapply(gtsetGenesList, function(gr) { unlist(strsplit(as.character(gr), ",", fixed=TRUE))})
	gtsetGenesList <- lapply(gtsetGenesList, function(gr) { 
						 	tmp <- unlist(strsplit(as.character(gr), ",", fixed=TRUE))
						 	gsub(pattern=" ", replacement="",  x=tmp, fixed=TRUE) # Remove spaces (4David).
						 })
	 
	names(gtsetGenesList) <- paste(as.character(geneTermSets[,1]), ".", unlist(lapply(table(factor(geneTermSets[,1]))[unique(geneTermSets[,1])], function(x) 1:x)), sep="")
	 
	nGenes <- length(unique(unlist(gtsetGenesList)))

	# Genes - GeneTermSet matrix
	gtSetGenesMatrix <- matrix(ncol=nGtSets, nrow=nGenes, data=0)
	colnames(gtSetGenesMatrix) <- names(gtsetGenesList)
	rownames(gtSetGenesMatrix) <- sort(unique(unlist(gtsetGenesList)))
	for( gts in  names(gtsetGenesList) )
	{
		gtSetGenesMatrix[eval(gtsetGenesList[[gts]]),gts]<-1
	}

	# Metagroups / Clusters 
	clusterGenesList <- list()
	for(cl in unique(geneTermSets[,1]))
	{
	 clusterGenesList[[cl]] <- geneTermSets[which(geneTermSets[,1] == cl),"Genes"]
	 clusterGenesList[[cl]] <- unlist(lapply(clusterGenesList[[cl]], function(gr) { strsplit(as.character(gr), ",", fixed=TRUE)})) 
	 clusterGenesList[[cl]] <- gsub(pattern=" ", replacement="",  x=clusterGenesList[[cl]], fixed=TRUE)
	}
	
	nGroups <- length(clusterGenesList)
	nGenes <- length(unique(unlist(clusterGenesList)))
	metagroupGenesMatrix <- matrix(ncol=nGroups, nrow=nGenes, data=0)
	colnames(metagroupGenesMatrix) <- names(clusterGenesList)
	rownames(metagroupGenesMatrix) <- sort(unique(unlist(clusterGenesList)))
	for( mg in  names(clusterGenesList) )
	{
		metagroupGenesMatrix[eval(clusterGenesList[[mg]]),mg]<-1
	}
	#colnames(metagroupGenesMatrix) <- paste(grType, names(clusterGenesList), sep="")

	return(list(metagroupGenesMatrix=metagroupGenesMatrix, gtSetGenesMatrix=gtSetGenesMatrix, filteredOut=filtrar))
}
