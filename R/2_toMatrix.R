
#tables <- toTables_gtLinker( results$globalMetagroups, results$geneTermSets)
#tables <- toMatrix(results$geneTermSets, results$globalMetagroups[,"Silhouette Width"])

#attribute: Sorted by metagroup/Cluster!!

toMatrix <- adjMatrix <- function(results, attribute=NULL, threshold=0, key="Genes", removeFiltered=NULL)
{
	# Check arguments
	grType <- "Mg"
	metagroups <- NULL
	if(is.list(results) && ("clusters" %in% names(results))) grType <- "Cl"
	if(is.list(results) && ("metagroups" %in% names(results))) metagroups <- results$metagroups   # Para david no hace falta, son lo mismo
	if(is.list(results) && ("geneTermSets" %in% names(results))) geneTermSets <- results$geneTermSets
	if(is.data.frame(results)) geneTermSets <- results
	if(!is.data.frame(geneTermSets)) stop("results should be list returned by getResults_gtLinker() or getResults_David().")
	if(dim(geneTermSets)[1] == 0) stop("0 gene term sets.")
	if(!is.numeric(threshold)) stop("threshold should be a number.")
	if(!is.null(attribute) && !is.data.frame(attribute)) stop("attribute should be the data.frame returned by getResults...().")	
		if(is.null(attribute) && threshold!=0) stop("To filter provide an attribute.")
	key <- capitalize(key) 
	if(!key %in% c("Genes", "Terms")) stop("key should be either 'Genes' or 'Terms'.")
	if(is.null(removeFiltered)) removeFiltered <- TRUE
	if(!is.logical(removeFiltered)) stop("removeFiltered should be either TRUE or FALSE.")
	if(removeFiltered && is.null(metagroups))  removeFiltered <- FALSE
	if(removeFiltered && key!="Terms") removeFiltered <- FALSE
	
	# Initialize
	if(key == "Genes") sepChar <- ","
	if(key == "Terms") sepChar <- ";"
	
	cols <- which(colnames(geneTermSets) %in% c("Cluster", "Metagroup", key))
	geneTermSets <- as.matrix(geneTermSets[,cols, drop=FALSE]) # Group & Genes
	
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
		
		# Filter MG:
		if(!is.null(metagroups)) metagroups <- metagroups[which(!rownames(metagroups) %in% filtrar),]
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
	gtSetsList <- as.character(geneTermSets[,key])
	gtSetsList <- lapply(gtSetsList, function(gr) { 
						 	unlist(strsplit(as.character(gr), sepChar, fixed=TRUE))
						 })
	 
	# Columna 1: Metagrupo/Cluster al que pertenece
	names(gtSetsList) <- paste(as.character(geneTermSets[,1]), ".", unlist(lapply(table(factor(geneTermSets[,1]))[unique(geneTermSets[,1])], function(x) 1:x)), sep="")
	 
	# Remove filtered terms from metagroups (only GeneTermLinker)
	if(removeFiltered)
	{
		nonFilteredList <- as.character(metagroups[,key])
		names(nonFilteredList) <- rownames(metagroups)
		nonFilteredList <- lapply(nonFilteredList, function(gr) { 
			unlist(strsplit(as.character(gr), sepChar, fixed=TRUE))
		})
		
		newGtSetsList <- NULL
		for(mg in names(nonFilteredList))
		{
			gtsets <- grep(paste("^",mg, "[.]", sep=""), names(gtSetsList))
			gtSetsList[gtsets] <-	lapply(gtSetsList[gtsets], function(gtset) gtset[gtset %in% nonFilteredList[[mg]]])
			tmp <- unique(gtSetsList[gtsets])
			names(tmp) <- paste(mg, ".", 1:length(tmp), sep="")
			newGtSetsList <- c(newGtSetsList, tmp)
		}
		gtSetsList <- newGtSetsList
	}
	
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
	for(mg in unique(geneTermSets[,1])) # Non-filtered mg
	{
		mgList[[mg]] <- unlist(gtSetsList[grep(paste("^",mg, "[.]", sep=""), names(gtSetsList))])
		mgList[[mg]] <- unique(unlist(lapply(mgList[[mg]], function(mg) { strsplit(as.character(mg), sepChar, fixed=TRUE)}))) 
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
