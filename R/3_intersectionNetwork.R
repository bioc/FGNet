
intersectionNetwork <- function(metagroupsMatrix, plotType="dynamic", vLayout="kk", returnGraph=FALSE, vSize=12, vLabelCex=2/3, legendMg=NULL, grPrefix="", plotTitle="Nodes in several metagroups", keepColors=TRUE, plotAllMg=FALSE)
{
	#################################### Check arguments ################################################
	filteredOut <- NULL
	if(is.list(metagroupsMatrix) && (any(names(metagroupsMatrix) %in% c("metagroupsMatrix", "clustersMatrix"))))
	{
		if("filteredOut" %in% names(metagroupsMatrix)) filteredOut <- metagroupsMatrix$filteredOut
		
		if("metagroupsMatrix" %in% names(metagroupsMatrix)) metagroupsMatrix <- metagroupsMatrix$metagroupsMatrix
		if("clustersMatrix" %in% names(metagroupsMatrix)) metagroupsMatrix <- metagroupsMatrix$clustersMatrix		
	}
	if(!is.matrix(metagroupsMatrix))  stop("metagroupsMatrix should be the result returned by toMatrix().")
	
	if(!is.character(plotType))  stop('plotType should be either "static", "dynamic" or "none".') 
	plotType <- tolower(plotType)
	if(!plotType %in% c("none", "static", "dynamic", "withpng")) stop('plotType should be "static" or "dynamic".') 
	
	if(!is.logical(returnGraph)) stop("returnGraph should be TRUE or FALSE.")
	
	if(!is.numeric(vSize))  stop("vSize should be numeric.")
	if(!is.numeric(vLabelCex))  stop("vLabelCex should be numeric.")
	if(!is.character(plotTitle)) stop("plotTitle should be character.")
	if(!is.character(grPrefix)) stop("grPrefix should be character.")
	
	
	# vLayout can be either c("kamada", "circle","sugiyama") or a set layout
	if(any(!vLayout %in% c("kk", "circle", "sugiyama"))) 
	{
		if(!is.matrix(vLayout) || (dim(vLayout)[2] != 2)) stop('vLayout should be "sugiyama", "kk" or "circle".')
		vLayout <- list(vLayout)
	} 
	
	if(!is.null(legendMg))
	{
		if(!is.character(legendMg)) stop("legendMg should be the metagroups legend (i.e. main function).")
		if(dim(metagroupsMatrix)[2] != length(legendMg)) stop("The number of metagroups in the matrix and the legend do not match.")
	}
	
	# Initialize
	nMetagroups <- dim(metagroupsMatrix)[2]	
	
	######## Set colors
	if(keepColors==TRUE)
	{
		colores <- setColors(as.character(sort(as.numeric(c(colnames(metagroupsMatrix), filteredOut)))))[colnames(metagroupsMatrix)]
	}else
	{
		colores <- setColors(colnames(metagroupsMatrix))
	}
	
	# Substitute terms for they description (for the graph)
	if(grepl(":", rownames(metagroupsMatrix)[1]))
	{
		clMgMxNames <- clDescriptions(rownames(metagroupsMatrix))
		rownames(metagroupsMatrix) <- paste(clMgMxNames[,"Description"], " (", clMgMxNames[,"Annotation"], ")", sep="")
	}
	
	#####################################################################################################
	#################################### Create Matrices  ###############################################	
	genesInManyMg <- which(apply(metagroupsMatrix, 1, sum)>1)
	
	if(!plotAllMg) metagroupsMatrix <- metagroupsMatrix[names(genesInManyMg),apply(metagroupsMatrix[names(genesInManyMg),], 2, sum)>0]

	mgJoinedGraph <- NULL
	if((length(genesInManyMg)>1) && any(plotType %in% c("static", "dynamic")))
	{
		mgJoined <- metagroupsMatrix[names(genesInManyMg),]
		
		# Replace 1 by mg number (for vertex color)
		for(col in 1:ncol(mgJoined))
		{
			thisCol <- which(mgJoined[,col]==1)
			mgJoined[thisCol,col] <- colnames(mgJoined)[col]
		}	
		
		if(!is.null(legendMg))
		{
			colnames(mgJoined) <- paste(colnames(mgJoined), legendMg, sep="\n")
		}
		if(grPrefix!="") colnames(mgJoined) <- paste(grPrefix, colnames(mgJoined), sep=" ")
		
		# Create adjacency matrix (add nodes and mg to both, columns and rows)
		mgXmg <- matrix(data=0, ncol=dim(mgJoined)[2], nrow=dim(mgJoined)[2])
		rownames(mgXmg) <- colnames(mgJoined)
		
		genXgen <- matrix(data=0, ncol=dim(mgJoined)[1], nrow=dim(mgJoined)[1]+dim(mgJoined)[2])
		colnames(genXgen) <- rownames(mgJoined)
		
		mgJoined <- rbind(mgJoined, mgXmg)
		mgJoined <- cbind(genXgen, mgJoined)
		
		# Vertex colors (genes: white, mg:their color)
		vCols2 <- c(rep("white", length(genesInManyMg)), colores[colnames(metagroupsMatrix)])
	
		# Create graph
		mgJoinedGraph <- graph.adjacency(mgJoined, weighted=TRUE,   mode="undirected", diag=FALSE) 
		
		for(layoutName in vLayout)
		{
			if(!is.matrix(layoutName)) 
			{
				if(layoutName %in% "kk")
				{
					vertexLayout <- layout.kamada.kawai(mgJoinedGraph)
					layoutName <- "Kamada Kawai"
				}
				if(layoutName %in% "circle")
				{
					vertexLayout <- layout.circle(mgJoinedGraph)
					layoutName <- "Circle"
				}
				if(layoutName %in% "sugiyama") 
				{ 
					vertexLayout <- layout.sugiyama(mgJoinedGraph)$layout
					layoutName <- "Sugiyama"
				}
			}else # Else: Layout is already provided
			{
				vertexLayout <- layoutName
				layoutName <- "Personalized"
			}
			
			if(plotType =="dynamic")	
			{
				tkplot(mgJoinedGraph, layout=vertexLayout, vertex.color=vCols2, edge.color=colores[as.character(E(mgJoinedGraph)$weight)], edge.width=2, vertex.size=vSize, vertex.label.cex=vLabelCex)
			}else
			{
				if(plotType!= "none")
				{
					plot(mgJoinedGraph, layout=vertexLayout, vertex.color=vCols2, edge.color=colores[as.character(E(mgJoinedGraph)$weight)], edge.width=2, vertex.size=vSize, vertex.label.cex=vLabelCex)
					title(main=plotTitle, sub=paste("Layout: ", layoutName, sep=""))
				}
			}
		}
	}else{
		if(length(genesInManyMg)<=1) warning("There is nothing to plot. There is no intersection between metagroups/clusters.")
	}
	
	if(returnGraph) return(mgJoinedGraph)
}
