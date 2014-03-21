

# Quitar: metagroupAttribute=NULL --> columnas en metagroupsMatrix: ordenadas por el atributo
functionalNetwork <- function(incidMatrices, plotType="static", returnGraph=FALSE, plotTitle="Functional Network", vSize=12, vLabelCex=3/4, vLayout=NULL, bgTransparency=0.4, legendMg=NULL, keepColors=TRUE, eColor="#323232", weighted=FALSE)
{
	#####################################################################################################
	#################################### Check arguments ################################################
	filteredOut <- NULL
	gtSetsMatrix <- NULL
	
	grType <- "Mg"
	if(is.list(incidMatrices) && (names(incidMatrices)[1] == "clustersMatrix"))   grType <- "Cl"
	
	if(is.matrix(incidMatrices)) metagroupsMatrix <- incidMatrices
	if(is.list(incidMatrices) && (names(incidMatrices) %in% c("metagroupsMatrix", "clustersMatrix")))
	{
		if(names(incidMatrices)[1] == "clustersMatrix") metagroupsMatrix <- incidMatrices$clustersMatrix
		if(names(incidMatrices)[1] == "metagroupsMatrix") metagroupsMatrix <- incidMatrices$metagroupsMatrix
		if("gtSetsMatrix" %in% names(incidMatrices))  gtSetsMatrix <- incidMatrices$gtSetsMatrix
		if("filteredOut" %in% names(incidMatrices))   filteredOut <- incidMatrices$filteredOut
	}
		
	if(is.matrix(metagroupsMatrix) && is.null(gtSetsMatrix))
	{
		gtSetsMatrix <- metagroupsMatrix
		filteredOut <- NULL
		
		warning("Only metagroupsMatrix provided.")
		if(plotTitle=="Functional Network") plotTitle <- "Metagroups Network"
	}	
	
	
	if(!is.matrix(metagroupsMatrix))  stop("metagroupsMatrix should be the result returned by toMatrix().")
	if(!is.matrix(gtSetsMatrix))  stop("gtSetsMatrix should be the result returned by toMatrix().")
	
	if(!is.character(plotType))  stop('plotType should be either "static", "dynamic" or "none".') 
	plotType <- tolower(plotType)
	if(!plotType %in% c("static", "dynamic", "withpng", "none")) stop('plotType should be in "static", "dynamic" or "none".') 
	
	if(!is.logical(returnGraph)) stop("returnGraph should be TRUE or FALSE.")
	
	if(!is.numeric(vSize))  stop("vSize should be numeric.")
	if(!is.numeric(vLabelCex))  stop("vLabelCex should be numeric.")
	if(!is.null(bgTransparency) && !is.numeric(bgTransparency))  stop("bgTransparency should be numeric.")
	if(!is.character(plotTitle)) stop("plotTitle should be character.")
	
	if(!is.null(vLayout) && (!is.matrix(vLayout) || (dim(vLayout)[2] != 2))) stop("Non-valid layout.")
	# eColor
	if(!is.logical(weighted)) stop("weighted should be either TRUE or FALSE.")
	
	plotLegend <- TRUE
	if(is.logical(legendMg) && legendMg==FALSE)
	{
		plotLegend <- FALSE
		legendMg <- NULL
	}
	if(!is.null(legendMg))
	{
		if(!is.character(legendMg)) stop("legendMg should be the metagroups legend (i.e. main function).")
		if(dim(metagroupsMatrix)[2] != length(legendMg)) stop("The number of metagroups in the matrix and the legend do not match.")
		legendMg <- paste(":", legendMg, sep="")
	}	
	
	
	# Initialize
	nGenes <- dim(metagroupsMatrix)[1]
	nMetagroups <- dim(metagroupsMatrix)[2]
	
	# Substitute terms for they description (for the graph)
	if(grepl(":", rownames(metagroupsMatrix)[1]))
	{
		clMgMxNames <- clDescriptions(rownames(metagroupsMatrix))
		rownames(metagroupsMatrix) <- paste(clMgMxNames[,"Description"], " (", clMgMxNames[,"Annotation"], ")", sep="")
		
		clGtSetMxNames <- clDescriptions(rownames(gtSetsMatrix))
		rownames(gtSetsMatrix) <- paste(clGtSetMxNames[,"Description"], " (", clGtSetMxNames[,"Annotation"], ")", sep="")
	}
	
	#####################################################################################################
	####################################  Plot settings  ################################################
	
	######## Set Metagroup colors
	if(keepColors==TRUE)
	{
		if(all(!is.na(as.numeric(colnames(metagroupsMatrix)))))
		{
			allMgSorted <- as.character(sort(as.numeric(c(colnames(metagroupsMatrix), filteredOut))))
		} else {
			allMgSorted <- sort(as.character(c(colnames(metagroupsMatrix), filteredOut)))
		}
		
		colores <- setColors(allMgSorted)[colnames(metagroupsMatrix)]
		trCols <- setColors(allMgSorted, transparency=bgTransparency)[colnames(metagroupsMatrix)]
	}else
	{
		colores <- setColors(colnames(metagroupsMatrix))
		trCols <- setColors(colnames(metagroupsMatrix), transparency=bgTransparency)
	}
	
	# Set colors for Nodes in only a group 
	vColor <- apply(metagroupsMatrix, 1, function(x) {
		if (sum(x) != 1)  return("#FFFFFF")
		else
		{
			return(colores[which(x ==1)])
		}
	})
	# Metagroup background color
	markGroup <- apply(metagroupsMatrix, 2, function(x) which(x==1))	
	if(!is.list(markGroup)) markGroup <- split(markGroup, rep(1:ncol(markGroup), each = nrow(markGroup))) #in case of same number of vertex in each group...
	
	
	#####################################################################################################
	#################################### Create Matrices  ###############################################
	# Gene - Gene edges 
	nCommonGTsets <- matrix(ncol=nGenes, nrow=nGenes, data=0)
	rownames(nCommonGTsets) <- rownames(metagroupsMatrix)
	colnames(nCommonGTsets) <- rownames(metagroupsMatrix)
	# Counts
	nCommonMgroups<- nCommonGTsets
	# Fill:
	for(gen1 in rownames(metagroupsMatrix))
		for(gen2 in rownames(metagroupsMatrix))
		{
			if(gen1!=gen2) nCommonMgroups[gen1, gen2] <- sum(metagroupsMatrix[gen1,] + metagroupsMatrix[gen2,]==2)
			if(gen1!=gen2) nCommonGTsets[gen1, gen2] <- sum(gtSetsMatrix[gen1,] + gtSetsMatrix[gen2,]==2)
		}
	#if(gen1!=gen2) edgeMatrix[gen1, gen2] <- any(abs(metagroupsMatrix[gen1,] + metagroupsMatrix[gen2,])==2) + any(abs(gtSetsMatrix[gen1,] + gtSetsMatrix[gen2,])==2)	
	
	#####################################################################################################
	####################################   Create Graph   ###############################################
	# weighted=TRUE The name of the edge attribute will be weight.
	# weighted=NULL (Unweighted graph. The elements of the adjacency matrix give the number of edges between vertices.)
	adjCommonEdges <- graph.adjacency(nCommonGTsets, weighted=TRUE, mode="undirected", diag=FALSE) 
	
	######## Layout
	# weighted=NULL (Unweighted graph. The elements of the adjacency matrix give the number of edges between vertices.) -> No parece importar para el layout
	vertexLayout <- vLayout
	if(is.null(vertexLayout))
	{
		graph4layout <- graph.adjacency(nCommonMgroups, weighted=NULL,   mode="undirected", diag=FALSE) 
		if(all(vColor=="#FFFFFF"))
		{
			# Parecido a fruchterman reingold, pero apelotona un poco menos los nodos comunes
			vertexLayout <- layout.kamada.kawai(graph4layout)	
		}else
		{
			# Si son dos grupos separados (no tienen ningun nodo en comun) kamada kawai superpone los dos grupos = MALO
			vertexLayout <- layout.fruchterman.reingold(graph4layout)#,  area=200000, repulserad=200000*vcount(graph4layout))
		}
	}
	
	#####################################################################################################
	#########################################  Draw PLOTS  ##############################################
	
	if(weighted && ecount(adjCommonEdges)>0)
	{
		eWidth <- E(adjCommonEdges)$weight
		eWidth <- floor((eWidth/max(eWidth))*5)
		if(min(eWidth)<1) eWidth <- eWidth +1
		if(min(eWidth) == max(eWidth)) eWidth <- 1
	}else eWidth <- 1
		
	if(plotType == "withpng")
	{
		require(png)		
		png("nwPreview.png", width = 800, height = 800)
		plot(adjCommonEdges, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.frame.color="#555555", vertex.size=vSize, vertex.label.color="#000000", vertex.label.cex=vLabelCex, mark.groups=markGroup, mark.col=trCols, mark.border=colores)#,  ,  mark.expand=2, , mark.shape=1)			
		if(plotLegend) legend(-1.2, -1.2, legend=paste(grType, names(colores), legendMg, sep="")[order(as.numeric(names(colores)))], fill=colores[order(as.numeric(names(colores)))], bty="n", xjust=0, yjust=0)
		dev.off()
	}
	
	if(plotType =="dynamic")
	{
		if(ecount(adjCommonEdges)>0)
		{
			tkplot(adjCommonEdges, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.size=vSize, vertex.label.color="black", vertex.label.cex=vLabelCex, vertex.frame.color=sub("#FFFFFF", "#888888", vColor)) 
		} else
		{
			warning("The dynamic plot is not available for networks without edges, try the static plot.") # Tkplot error: Error in mapply(function(from, to, id) .tkplot.create.edge: zero-length inputs cannot be mixed with those of non-zero length
		}
		intersectionNetwork(metagroupsMatrix, plotType="dynamic", vLayout="kk", returnGraph=FALSE, vSize=vSize, vLabelCex=vLabelCex, legendMg=legendMg)
	}else
	{
		# PLOT
		if(plotType !="none")
		{
			plot(adjCommonEdges, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.frame.color="#555555", vertex.size=vSize, vertex.label.color="#000000", vertex.label.cex=vLabelCex, mark.groups=markGroup, mark.col=trCols, mark.border=colores)#,  ,  mark.expand=2, , mark.shape=1)		
			if(plotLegend) legend(-1.2, -1.2, legend=paste(grType, names(colores)[order(as.numeric(names(colores)))], legendMg, sep=""), fill=colores[order(as.numeric(names(colores)))], bty="n", xjust=0, yjust=0)
			title(plotTitle)
		}
	}
	
	if(returnGraph) return(adjCommonEdges)
}


