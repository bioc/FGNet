
# Returns a list with the metagroup term names and links
getMGTerms <- function(globalMetagroups, grType, org=NULL)
{
	if(is.null(org)) org <- "map"	
	if(tolower(org)=="sc") org <- "sce"	
	if(tolower(org)=="hs") org <- "hsa"	
	
	mgTerms<- sapply(globalMetagroups$Terms, function(x) strsplit(as.character(x), split=";"))
	names(mgTerms)<-paste(grType,rownames(globalMetagroups), sep=" ")
	
	termsDescriptions <- list()		
	goIds <- list()
	for(mgName in names(mgTerms))
	{
		descripciones <- clDescriptions(mgTerms[[mgName]], org)
		
		descripciones <- descripciones[order(descripciones[,"Description"]),, drop=FALSE]
		termsDescriptions[[mgName]] <- descripciones
		
		# GO
		goIds[[mgName]] <- unlist(apply(descripciones,1, function(x) if(x["Annotation"] =="GO") return(x["goID"])))	# GO ID is saved in fourth slot 
	}
	return(list(termsDescriptions=termsDescriptions, goIds=goIds))
}

#representativeTerm(mgTerms$termsDescriptions)

representativeTerm <- function(termsDescriptions)
{
	
	# Search for a representative term (based only on key words! not function!)
	reprTerms <- rep(NA, length(termsDescriptions))
	for(index in 1:length(termsDescriptions))
	{
		termsMg <- as.matrix(lapply(termsDescriptions, function(x) cbind(x[,"Description"]))[[index]])
		
		rownames(termsMg) <- NULL
		splittedTerms <- unlist(strsplit(termsMg, " "))
		splittedTerms <- unlist(strsplit(splittedTerms, "-", fixed=TRUE))
		splittedTerms <- unlist(strsplit(splittedTerms, "/", fixed=TRUE))
		splittedTerms <- unlist(strsplit(splittedTerms, ",", fixed=TRUE))
		splittedTerms <- unlist(strsplit(splittedTerms, ".", fixed=TRUE))
		splittedTerms <- unlist(strsplit(splittedTerms, "(", fixed=TRUE))
		splittedTerms <- unlist(strsplit(splittedTerms, ")", fixed=TRUE))
		
		sortedFreq <- sort(table(tolower(splittedTerms)), decreasing=TRUE)
		
		commonWords <- sortedFreq[sortedFreq >= quantile(sortedFreq, 0.9)]
		commonWords <- commonWords[grepl("[[:alnum:]]", names(commonWords))]
		
		selectedTerm <- NULL
		if(any(tolower(termsMg) %in% names(commonWords)))
		{
			selectedTerm <- names(commonWords)[which(names(commonWords) %in% tolower(termsMg))[1]] # Select the one with highest freq.
			
			if(commonWords[1] >= commonWords[selectedTerm]*2) selectedTerm <- NULL 	
			selectedTerm <- capitalize(selectedTerm)
		}
		if(is.null(selectedTerm))
		{
			termsWithCommonWords <- sapply(names(commonWords), function(x) grep(x, tolower(termsMg)))
			countTermsWords <- table(unlist(termsWithCommonWords)) #Number of common words per term		
			selectedTerms <- termsMg[as.numeric(names(countTermsWords[countTermsWords == max(countTermsWords)])),] # Terms with most common words
			ncharTerms <- nchar(selectedTerms)
			selectedTerm <- selectedTerms[ncharTerms==min(ncharTerms)][1]
		}
		
		reprTerms[index] <- selectedTerm
		if(nchar(reprTerms[index])>30) reprTerms[index] <- paste(substr(reprTerms[index], 1, 30), "-", sep="")
	}
	return(reprTerms)
}

# (Create html links)
buildTermsTable <- function(termsDescriptions)
{
	termLinks <- lapply(termsDescriptions, function(mgMx) {
		if(!is.matrix(mgMx)) mgMx <- t(as.matrix(mgMx)) # 1 row?
		descrLinks <- apply(mgMx, 1, function(term){
			if(is.na(term["Link"]))
			{
				link <- term["Description"]
			}else
			{
				link <- paste("<a href='", term["Link"],"' target='_blank'>", term["Description"] ,"</a>", sep="")
			}
			return (link)
		})
		
		descrLinks <- cbind(descrLinks, mgMx[,"Annotation"])
		colnames(descrLinks) <- NULL
		rownames(descrLinks) <- NULL
		return(descrLinks)
	})
	
	return(termLinks)
}



# Returns the link to the ontology tree for each metagroup go terms
goTreeLinks <- function(goIds, folder, downloadGOtree)
{
	goTerms <- sapply(goIds, function(x) if(length(x)>0) paste("%22GO%3A", x,"%22%3A{%22fill%22%3A%22%23ccccff%22}", sep="",collapse=","))
	if(length(goTerms)>0)
	{
		goLinks <- paste("http://amigo.geneontology.org/visualize?term_data={",goTerms ,"}&term_data_type=json&mode=amigo&format=png&inline=false", sep="")
		names(goLinks) <- names(goIds)
		
		fileNames <- paste(folder, paste("GO_", gsub("s ","_",names(goIds)),".png", sep=""),sep="")
		names(fileNames) <- names(goIds)
		
		if(downloadGOtree)
		{
			for(i in 1:length(goLinks)) 
			{
				fileNames[i] <- tryCatch({
						download.file(url=goLinks[i], destfile=fileNames[i], quiet=TRUE)
						fileNames[i] # se ha descargado bien: linkar
				}, error = function(e) {
						if (file.exists(fileNames[i])) file.remove(fileNames[i])
					  return(goLinks[i]) # No se ha podido descargar. Link a la URL
				})
			}
		}else{
			fileNames <- goLinks
		}
		return(fileNames)
	} else 
	{
		return("")
	}
}


# Main function
createHtml <- function(htmlFileName, results, tablesGenes, tablesTerms, metagroupAttributeName, threshold, jobID=NULL, organism=NULL, annotations=NULL, genes=NULL, serverWeb, downloadGOtree=TRUE)
{
	#####################################################################################################
	####################################   Initializations   ############################################
	
	tables <- tablesGenes
	grType <- names(results)[1] # $metagroups or $clusters
	
	# David or gtLinker?
	if(grType == "metagroups") 
	{
		grType <- "Metagroup"
		grPrefix <- "mg"
		
		mgIncMatrix <- tables$metagroupsMatrix
		termsNumMg <- apply(tablesTerms$metagroupsMatrix,1,sum)
	}
	if(grType == "clusters") 	
	{
		grType <- "Cluster"
		grPrefix <- "cl"
		
		mgIncMatrix <- tables$clustersMatrix
		termsNumMg <- apply(tablesTerms$clustersMatrix,1,sum)
	}
	
	rawMetagroups <- results[[1]]
	filteredOut <- tables$filteredOut
	
	nGenesTotal <- length(unique(unlist(sapply(rawMetagroups[,"Genes"], function(x) strsplit(as.character(x), split=",")))))
	nGenesInNotFilteredMg <- dim(tables[[1]])[1]
	nRawMg <- dim(rawMetagroups)[1]
	
	# Get metagroup colors
	colores <- setColors(as.character(sort(as.numeric(c(colnames(mgIncMatrix), filteredOut)))))[colnames(mgIncMatrix)]
	
	globalMetagroups <- rawMetagroups[colnames(mgIncMatrix),]
	
	#### Set locations
	# Folter to save images etc: Same as downloaded results	
	rawResults <- results$fileName # Whole path
	tmp <- gregexpr("./.", rawResults)[[1]]
	folder <- substring(rawResults, first=1, last=tmp[length(tmp)]+1) # Location (Folder)
	
	if(folder!="") rawResults <- substring(rawResults, first=attr(gregexpr(folder, rawResults)[[1]], "match.length")+1) # Only file name
	
	tmp <- gregexpr("/", htmlFileName)[[1]]
	htmlRoot <- substring(htmlFileName, first=1, last=tmp[length(tmp)])
	
	# Is the folder a subfolder of the html root? (-> make HTML routes relative)
	folderInRoot <- gregexpr(htmlRoot, folder)[[1]]
	if(folderInRoot[1] > 0)
	{
		folder <- substring(folder, first=attr(folderInRoot, "match.length")+1)
	}	
	
	# Get metagroup terms
	mgTerms <- getMGTerms(globalMetagroups, grType, organism)
	
	# Add asterisc to terms in several metagroups:
	termsMultipleMg <- termsNumMg[which(termsNumMg>1)]
	if(length(termsMultipleMg)>0)
	{
		mgTerms$termsDescriptions <- sapply(mgTerms$termsDescriptions, function(x) 
		{
			tmp <- rownames(x) %in% names(termsMultipleMg)
			if(any(tmp))  x[tmp,"Description"] <- paste(x[tmp,"Description"], " *", sep="")
			x
		})
	}
		
	# Terms table to write in html: 
	termsTables <- buildTermsTable(mgTerms$termsDescriptions)
	
	goIds <- mgTerms$goIds
	goLinks <- goTreeLinks(goIds,folder, downloadGOtree)																					
	
	# Copiar CSS a la carpeta actual...
	cssDir <- file.path(system.file('css', package='FGNet'))
	cssFile <- paste(cssDir, "functionalNetworks.css", sep="/")
	if(folder=="") file.copy(cssFile, "./")
	if(folder!="") file.copy(cssFile, folder)
	
	# Add representative term to DAVID legend
	legendMg <- NULL
	if(grType ==  "Cluster") legendMg <- paste(representativeTerm(mgTerms$termsDescriptions), ", ...", sep="")
	
	#####################################################################################################
	#################################### Create Plots  ###################################################
	
	networkPlot <- paste(folder, "nwFunctionalNetwork.png", sep="")
	networkPlot2a <- paste(folder, "nwIntersection_kk.png", sep="")
	networkPlot2b <- paste(folder, "nwIntersection_circle.png", sep="")
	#networkPlot2c <- paste(folder, "nwIntersection_sugiyama.png", sep="")	
	networkTerms <- paste(folder, "nwTerms.png", sep="")
	distancePlot <- paste(folder, "plot_Distance.png", sep="")
	iGraphFile <- paste(folder, "iGraph.RData", sep="")
	simplifiedTable <- paste(folder, "simplifiedTermsTable.txt", sep="")
	
	# Set node and label size
	numNodes <- sum(dim(mgIncMatrix))
	if(numNodes < 150)
	{
		vSize <- 12
		vLabelCex <- 0.9
	}
	if(numNodes >= 150)
	{
		vSize <- 10
		vLabelCex <- 3/4
	}	
	if(numNodes > 300)
	{
		vSize <- 8
		vLabelCex <- 2/3
	}	
	
	png(networkPlot, width = 800, height = 800)
	iGraph <- functionalNetwork(tables, plotType="static", vSize=vSize, vLabelCex=vLabelCex, returnGraph=TRUE, legendMg=legendMg)
	dev.off()
	png(networkPlot2a, width = 800, height = 800)
	intersectionGraph <- intersectionNetwork(tables, vLayout=c("kk"), plotType="static", vSize=vSize, vLabelCex=vLabelCex, grPrefix=grPrefix, returnGraph=TRUE)
	dev.off()
	if(!is.null(intersectionGraph))
	{
		png(networkPlot2b, width = 800, height = 800)
		intersectionNetwork(tables, vLayout=c("circle"), plotType="static",  vSize=vSize, vLabelCex=vLabelCex, grPrefix=grPrefix)
		dev.off()
# 		png(networkPlot2c, width = 800, height = 800)
# 		intersectionNetwork(tables, vLayout=c("sugiyama"), plotType="static",  vSize=vSize, vLabelCex=vLabelCex, grPrefix=grPrefix)
# 		dev.off()
	}
	
	if(length(termsMultipleMg)>0)
	{
		png(networkTerms, width = 800, height = 800)
		intersectionNetwork(tablesTerms, plotType="static", plotTitle=paste("Terms in several ", tolower(grType), "s", sep=""),  plotAllMg=FALSE, vSize=vSize, vLabelCex=vLabelCex, grPrefix=grPrefix) 
		dev.off()
	}
	
	png(distancePlot, width = 600, height = 600)
	distanceMatrix <- plotMetagroupsDistance(tables)
	dev.off()
	
	save(iGraph, file=iGraphFile)
	
	#####################################################################################################
	#################################### Create HTML  ###################################################
	p=openPage(htmlFileName, link.css=paste(folder, "functionalNetworks.css", sep=""))
	
	# Header
	hwrite('Functional Gene Network', p, heading=1)
	
	# Query parameters
	hwrite("Parameters of the query: ", p, heading=2)
	if(!is.null(organism))
	{
		hwrite("Organism: ", p, class='InfoLabel', br=FALSE)
		if(decapitalize(organism)=="sc") organism <- "Sc (Saccharomyces cerevisiae)"
		if(decapitalize(organism)=="hs") organism <- "Hs (Homo Sapiens)"
		hwrite(organism, p, br=TRUE)
	}
	if(!is.null(annotations))
	{
		hwrite("Annotations: ", p, class='InfoLabel', br=FALSE)
		hwrite(paste(annotations, sep="", collapse=", "), p, br=TRUE)
		#hwrite(annotations, p, dim=c(length(annotations),1), border=0, br=FALSE)
	}
	if(!is.null(genes))
	{
		hwrite("Genes: ", p, class='InfoLabel', br=FALSE)
		hwrite(paste("(", length(genes), ") ", sep=""), p, br=FALSE)
		hwrite(paste(sort(genes), sep="", collapse=", "), p, br=TRUE)
	}
	hwrite("Server/Tool: ", p, class='InfoLabel', br=FALSE)
	hwrite(serverWeb, link=paste(serverWeb, '" target="_blank', sep=""),p, br=TRUE)
	hwrite("Raw results from functional enrichment and clustering (.txt): ", p,  br=FALSE)
	
	if(grType == "Metagroup") 
	{
		hwrite('[Global overview] ', p, link=paste(folder, rawResults, sep="" ), br=FALSE)
		hwrite(paste('[Mg', 1:nRawMg,']', sep=""), p, link=paste(folder, substring(rawResults, first=1, last=gregexpr("_global_overview.txt", rawResults)[[1]][1]), tolower(grType), "_", 1:nRawMg, ".txt", sep=""), table=FALSE, br=FALSE)	
		hwrite('',p, br=TRUE)
		if(!is.null(jobID))
		{
			hwrite("Job ID: ", p, class='InfoLabel', br=FALSE)
			hwrite(jobID, p, br=TRUE)
		}
	}
	if(grType == "Cluster") 
	{
		hwrite('[DavidClustering] ', p, link=paste(folder, rawResults, sep="" ), br=FALSE)
	}
	
	# Results header				
	hwrite("Results: ", p, heading=2)
	if(grType == "Metagroup") 
	{
		hwrite("Number of metagroups: ", p, class='InfoLabel', br=FALSE)
		hwrite(nRawMg, p, br=TRUE)
		hwrite("Number of genes included in all metagroups: ", p, class='InfoLabel', br=FALSE)
		hwrite(nGenesTotal, p, br=TRUE)
		hwrite("Filtered metagropus ", p, class='InfoLabel', br=FALSE)
		hwrite(paste(" (", metagroupAttributeName, " < ", threshold, "): ", sep=""), p, br=FALSE)
		if(length(filteredOut)>0)  hwrite(paste("Mg", filteredOut, sep="", collapse=", "), p, br=TRUE) 
		else hwrite("None", p, br=TRUE) 
		hwrite("Number of genes included in non-filtered metagroups: ", p, class='InfoLabel', br=FALSE)
		hwrite(nGenesInNotFilteredMg, p, br=TRUE)
	}
	if(grType == "Cluster") 
	{
		hwrite("Number of clusters: ", p, class='InfoLabel', br=FALSE)
		hwrite(nRawMg, p, br=TRUE)
		hwrite("Number of genes included in all clusters: ", p, class='InfoLabel', br=FALSE)
		hwrite(nGenesTotal, p, br=TRUE)
		hwrite("Filtered clusters ", p, class='InfoLabel', br=FALSE)
		hwrite(paste(" (", metagroupAttributeName, " < ", threshold, "): ", sep=""), p, br=FALSE)
		if(length(filteredOut)>0)  hwrite(paste("Cl", filteredOut, sep="", collapse=", "), p, br=TRUE) 
		else hwrite("None", p, br=TRUE) 
		hwrite("Number of genes included in non-filtered clusters: ", p, class='InfoLabel', br=FALSE)
		hwrite(nGenesInNotFilteredMg, p, br=TRUE)
	}
	hwrite("Functional network in other formats: ", p)
	#if(!is.null(pdfName)) hwrite('[PDF] ', p, link=paste(pdfName, sep="" ), br=FALSE)
	hwrite('[iGraph] ', p, link=paste(iGraphFile, sep="" ), br=TRUE)
	
	hwrite("Simplified metagroup-terms table (shown below as legend): ", p)
	hwrite('[Legend] ', p, link=paste(simplifiedTable, sep="" ), br=FALSE)
	
	# Functional Network plot
	#hwrite("Functional Network: ", p, heading=3)
	subImages <- ""
	closeButton <- '<a href="#close" title="Close" class="close">X</a>'
	if(!is.null(intersectionGraph))
	{		
		subImages <- c(hwrite(paste('Genes in several ', tolower(grType),'s: ', sep=""), class='InfoLabel', br=TRUE),
									 hwrite('(Different layouts)'),
									 
									 hwriteImage(networkPlot2a, class='intersectionNw', link='#nwIntersection_kk', br=FALSE),
									 hwrite(paste('<div id="nwIntersection_kk" class="modalDialog"><div>',closeButton, '<img border="0" src="',networkPlot2a,'" width="100%"></div></div>', sep="")),
									 									 
									 hwriteImage(networkPlot2b, class='intersectionNw', link='#nwIntersection_circle', br=FALSE),
									 hwrite(paste('<div id="nwIntersection_circle" class="modalDialog"><div>',closeButton, '<img border="0" src="',networkPlot2b,'" width="100%"></div></div>', sep="")))

		
		if(length(termsMultipleMg)>0) 
		{
			subImages <- c(subImages,
										 hwrite(paste('Terms in several ', tolower(grType),'s: ', sep=""), class='InfoLabel', br=TRUE),
										 hwrite('(Marked with * in the table)'),
										 hwriteImage(networkTerms, class='intersectionNw', link='#nwTerms', br=FALSE),
										 hwrite(paste('<div id="nwTerms" class="modalDialog"><div>',closeButton, '<img border="0" src="',networkTerms,'" width="100%"></div></div>', sep="")))
		}
	
	}
	imageTable <- c(hwriteImage(networkPlot, class='network', link='#functionalNetwork', br=FALSE),
									hwrite(paste('<div id="functionalNetwork" class="modalDialog"><div>',closeButton, '<img border="0" src="',networkPlot,'" width="100%"></div></div>', sep="")),
									hwrite(subImages, border=0, dim=c(length(subImages),1), class="ImageTable"))
	hwrite(c(imageTable), p, border=0, class="ImageTable")
	
	if(grType == "Metagroup") hwrite("Metagroups (sorted by Silhouette): ", p, br=TRUE)
	if(grType == "Cluster") hwrite("Clusters (sorted by Score): ", p, br=TRUE)
	# Terms table
	termsTable <- NULL
	txtTermsTable <- NULL
	for(mg in 1:length(termsTables))
	{
		mgName <- paste("<b>",names(termsTables)[mg],"</b>")
		attrs <- NULL
		if(grType == "Metagroup")
		{
			attrs <- c(paste("Silhouette: ", round(globalMetagroups[mg, "Silhouette Width"], 2), sep=""),
								 paste("P-value: ", signif(globalMetagroups[mg, "pValue"],2),sep=""))		
		}		
		if(grType == "Cluster")
		{	
			attrs <- c("", # "Silhouette  only in GtLinker
								 paste("Score: ", signif(as.numeric(as.character(globalMetagroups[mg, "EnrichmentScore"])),2),sep=""))
		}		
		tmpGenes <- paste(sort(unlist(strsplit(as.character(globalMetagroups[mg, "Genes"]), ","))), collapse=", ")
		linkGenes <- paste("Genes: ", '<a href="#close" title="', tmpGenes,'" class="tooltip"><span title="Genes">', globalMetagroups[mg, "nGenes"], '</span></a>',sep="")			 
		txtGenes <- paste(globalMetagroups[mg, "nGenes"], " genes: ", tmpGenes, sep="")
		attrsTxt <- c(attrs, txtGenes)
		attrs <- c(attrs, linkGenes)					 
		
		terms <- rbind(termsTables[[mg]], if(!is.na(goLinks[names(termsTables)[mg]])) cbind("", paste("<a href='",goLinks[names(termsTables)[mg]],"'' target='_blank' class='goLink'>[GO tree]</a>", sep="")))
		
		termsTable <- c(termsTable, hwrite(c(mgName, attrs), class='mgAttr', table=TRUE, border=0), 
										hwrite(terms, table=TRUE, col.class=list("termsRow", "annotRow")))
		txtTermsTable <- rbind(txtTermsTable, paste(c(names(termsTables)[mg], attrsTxt), collapse="\t"), rbind(mgTerms$termsDescriptions[[mg]])[, "Description", drop=FALSE])
		#cbind(mgTerms$termsDescriptions[[mg]][,"Description"]))
	}
	
	#bgCols <- c(strsplit(paste(substr(colores, 1, 7)[rownames(globalMetagroups)], "#FFFFFF", collapse=" "), split=" ")[[1]],"#FFFFFF")
	borderCols <- paste('border-color:', strsplit(paste(substr(colores, 1, 7), substr(colores, 1, 7), collapse=" "), split=" ")[[1]], sep="")
	hwrite(termsTable, p, border=1, class=rep(c('mgHeader', 'Terms'), length(termsTable)/2), row.style=borderCols, dim=c(length(termsTable),1))
	
	# Explanation of asterisc (there are terms in several metagroups):
	if(length(termsMultipleMg)>0) 
	{
		hwrite(paste("<i>* Terms marked with an asterisc are in several ", tolower(grType),"s.</i>", sep=""), p, br=TRUE)
		#hwriteImage(networkTerms, p, link=c(networkTerms), br=TRUE,border=0)	
	} else {
		hwrite(paste("<i>(There are no terms in more than one ", tolower(grType),")</i>", sep=""), p, br=TRUE)
	}
	
	# Distances Plot	
	if(!is.null(distanceMatrix))
	{
		hwrite(paste("Distances between ", grType,"s: ", sep=""), p, heading=3, br=TRUE)
		if(grType == "Metagroup")
		{
			colnames(distanceMatrix) <- paste("Mg", colnames(distanceMatrix), sep="")
			rownames(distanceMatrix) <- paste("Mg", rownames(distanceMatrix), sep="")
		}
		if(grType == "Cluster") 
		{
			colnames(distanceMatrix) <- paste("Cl", colnames(distanceMatrix), sep="")
			rownames(distanceMatrix) <- paste("Cl", rownames(distanceMatrix), sep="")
		}
		subImages <- c(hwrite('Distance matrix: ', class='InfoLabel', br=TRUE),
									 hwrite(round(distanceMatrix,2), table=TRUE, border=0, class='matrix', br=FALSE))
		imageTable <- c(hwriteImage(distancePlot, class='distance', link=distancePlot, br=TRUE), hwrite(subImages, border=0, dim=c(length(subImages),1)))
		
		hwrite(c(imageTable), p, table=TRUE, border=0) #, class="ImageTable"
	}
	closePage(p)
	
	write.table(txtTermsTable, file=simplifiedTable, row.names=FALSE, col.names=FALSE, quote=FALSE)
	
	browseURL(htmlFileName) # Open html
}
