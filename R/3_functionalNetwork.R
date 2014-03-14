
# metagroupGenesMatrix <- tables$metagroupGenesMatrix
# gtSetGenesMatrix <- tables$gtSetGenesMatrix

#functionalNetwork(tables$metagroupGenesMatrix, tables$gtSetGenesMatrix) 
#functionalNetwork(tables$metagroupGenesMatrix, tables$gtSetGenesMatrix, plotType="dynamic") 


# Quitar: metagroupAttribute=NULL --> columnas en MetagroupGenesMatrix: ordenadas por el atributo
functionalNetwork <- function(metagroupGenesMatrix, gtSetGenesMatrix=NULL, plotType="static", returnGraph=FALSE, plotTitle="Functional Network", vSize=12, vLabelCex=3/4, vLayout=NULL, bgTransparency=0.4, legendMg=NULL, keepColors=TRUE)
{
	# Libraries
	# require(igraph)
	# if (!library(igraph, logical.return=TRUE)) stop("Library igraph is required to plot the networks.")
	
	#####################################################################################################
	#################################### Check arguments ################################################
	filteredOut <- NULL
	
	if(is.matrix(metagroupGenesMatrix) && is.null(gtSetGenesMatrix))
	{
		gtSetGenesMatrix <- metagroupGenesMatrix
		filteredOut <- NULL
		warning("Only metagroupGenesMatrix provided.")
	}
	
	if(is.list(metagroupGenesMatrix) && (all(names(metagroupGenesMatrix) %in% c("metagroupGenesMatrix", "gtSetGenesMatrix", "filteredOut"))))
	{
		gtSetGenesMatrix <- metagroupGenesMatrix$gtSetGenesMatrix
		filteredOut <- metagroupGenesMatrix$filteredOut
		metagroupGenesMatrix <- metagroupGenesMatrix$metagroupGenesMatrix
	}
	if(!is.matrix(metagroupGenesMatrix))  stop("metagroupGenesMatrix should be the result returned by toMatrix().")
	if(!is.matrix(gtSetGenesMatrix))  stop("gtSetGenesMatrix should be the result returned by toMatrix().")
	
	if(!is.character(plotType))  stop('plotType should be either "static", "dynamic" or "none".') 
	plotType <- tolower(plotType)
	if(!plotType %in% c("static", "dynamic", "withpng", "none")) stop('plotType should be in "static", "dynamic" or "none".') 
	
	if(!is.logical(returnGraph)) stop("returnGraph should be TRUE or FALSE.")
	
	if(!is.numeric(vSize))  stop("vSize should be numeric.")
	if(!is.numeric(vLabelCex))  stop("vLabelCex should be numeric.")
	if(!is.null(bgTransparency) && !is.numeric(bgTransparency))  stop("bgTransparency should be numeric.")
	if(!is.character(plotTitle)) stop("plotTitle should be character.")
	
	if(!is.null(vLayout) && (!is.matrix(vLayout) || (dim(vLayout)[2] != 2))) stop("Non-valid layout.")
	
	plotLegend <- TRUE
	if(is.logical(legendMg) && legendMg==FALSE)
	{
		plotLegend <- FALSE
		legendMg <- NULL
	}
	if(!is.null(legendMg))
	{
		if(!is.character(legendMg)) stop("legendMg should be the metagroups legend (i.e. main function).")
		if(dim(metagroupGenesMatrix)[2] != length(legendMg)) stop("The number of metagroups in the matrix and the legend do not match.")
		legendMg <- paste(":", legendMg, sep="")
	}	
	
	# Initialize
	nGenes <- dim(metagroupGenesMatrix)[1]
	nMetagroups <- dim(metagroupGenesMatrix)[2]
	
	#####################################################################################################
	#################################### Create Matrices  ###############################################
	# Gene - Gene edges 
	nCommonGTsets <- matrix(ncol=nGenes, nrow=nGenes, data=0)
	rownames(nCommonGTsets) <- rownames(metagroupGenesMatrix)
	colnames(nCommonGTsets) <- rownames(metagroupGenesMatrix)
	# Counts
	nCommonMgroups<- nCommonGTsets
	# Fill:
	for(gen1 in rownames(metagroupGenesMatrix))
		for(gen2 in rownames(metagroupGenesMatrix))
		{
			if(gen1!=gen2) nCommonMgroups[gen1, gen2] <- sum(metagroupGenesMatrix[gen1,] + metagroupGenesMatrix[gen2,]==2)
			if(gen1!=gen2) nCommonGTsets[gen1, gen2] <- sum(gtSetGenesMatrix[gen1,] + gtSetGenesMatrix[gen2,]==2)
		}
	#if(gen1!=gen2) edgeMatrix[gen1, gen2] <- any(abs(metagroupGenesMatrix[gen1,] + metagroupGenesMatrix[gen2,])==2) + any(abs(gtSetGenesMatrix[gen1,] + gtSetGenesMatrix[gen2,])==2)	
	
	#####################################################################################################
	####################################  Plot settings  ################################################
	
	# Create graph
	# weighted=TRUE The name of the edge attribute will be weight.
	adjCommonEdges <- graph.adjacency(nCommonGTsets, weighted=TRUE,   mode="undirected", diag=FALSE) 
	
	######## Set colors
	if(keepColors==TRUE)
	{
		if(all(!is.na(as.numeric(colnames(metagroupGenesMatrix)))))
		{
			allMgSorted <- as.character(sort(as.numeric(c(colnames(metagroupGenesMatrix), filteredOut))))
		} else {
			allMgSorted <- sort(as.character(c(colnames(metagroupGenesMatrix), filteredOut)))
		}
		
		colores <- setColors(allMgSorted)[colnames(metagroupGenesMatrix)]
		trCols <- setColors(allMgSorted, transparency=bgTransparency)[colnames(metagroupGenesMatrix)]
	}else
	{
		colores <- setColors(colnames(metagroupGenesMatrix))
		trCols <- setColors(colnames(metagroupGenesMatrix), transparency=bgTransparency)
	}
	
	# Nodes in only a group 
	vColor <- apply(metagroupGenesMatrix, 1, function(x) {
		if (sum(x) != 1)  return("#FFFFFF")
		else
		{
			return(colores[which(x ==1)])
		}
	})
	# Metagroup Color
	markGroup <- apply(metagroupGenesMatrix, 2, function(x) which(x==1))	
	if(!is.list(markGroup)) markGroup <- split(markGroup, rep(1:ncol(markGroup), each = nrow(markGroup))) #in case of same number of vertex in each group...
	
	# Edge color
	eColor <- "#323232"	
	
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
	
	if(plotType == "withpng")
	{
		require(png)		
		png("nwPreview.png", width = 800, height = 800)
		plot(adjCommonEdges, layout=vertexLayout, edge.width=1, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.frame.color="#555555", vertex.size=vSize, vertex.label.color="#000000", vertex.label.cex=vLabelCex, mark.groups=markGroup, mark.col=trCols, mark.border=colores)#,  ,  mark.expand=2, , mark.shape=1)			
		if(plotLegend) legend(-1.2, -1.2, legend=paste("Mg", names(colores), legendMg, sep="")[order(names(colores))], fill=colores[order(names(colores))], bty="n", xjust=0, yjust=0)
		dev.off()
	}
	
	if(plotType =="dynamic")
	{
		tkplot(adjCommonEdges, layout=vertexLayout, edge.width=1, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.size=vSize, vertex.label.color="black", vertex.label.cex=vLabelCex, vertex.frame.color=sub("#FFFFFF", "#888888", vColor)) #edge.width=(E(adjCommonEdges)$weight)
		intersectionNetwork(metagroupGenesMatrix, plotType="dynamic", vLayout="sugiyama", returnGraph=FALSE, vSize=vSize, vLabelCex=vLabelCex, legendMg=legendMg)
	}else
	{
		# PLOT
		if(plotType !="none")
		{
			plot(adjCommonEdges, layout=vertexLayout, edge.width=1, edge.color=eColor, vertex.label=V(adjCommonEdges)$name, vertex.color=vColor, vertex.frame.color="#555555", vertex.size=vSize, vertex.label.color="#000000", vertex.label.cex=vLabelCex, mark.groups=markGroup, mark.col=trCols, mark.border=colores)#,  ,  mark.expand=2, , mark.shape=1)		
			if(plotLegend) legend(-1.2, -1.2, legend=paste("Mg", names(colores)[order(names(colores))], legendMg, sep=""), fill=colores[order(names(colores))], bty="n", xjust=0, yjust=0)
			title(plotTitle)
		}
	}
	
	if(returnGraph) return(adjCommonEdges)
}


