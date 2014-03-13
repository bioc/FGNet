
intersectionNetwork <- function(metagroupGenesMatrix, plotType="dynamic", vLayout="kk", returnGraph=FALSE, vSize=12, vLabelCex=2/3, legendMg=NULL, grPrefix="", plotTitle="Nodes in several metagroups", keepColors=TRUE)
{
	# Libraries
	# if (!library(igraph, logical.return=TRUE)) stop("Library igraph is required to plot the networks.")

	#################################### Check arguments ################################################
	filteredOut <- NULL
	if(is.list(metagroupGenesMatrix) && (all(names(metagroupGenesMatrix %in% c("metagroupGenesMatrix", "gtSetGenesMatrix", "filteredOut")))))
	{
		gtSetGenesMatrix <- metagroupGenesMatrix$gtSetGenesMatrix
		filteredOut <- metagroupGenesMatrix$filteredOut
		metagroupGenesMatrix <- metagroupGenesMatrix$metagroupGenesMatrix
	}
	if(!is.matrix(metagroupGenesMatrix))  stop("metagroupGenesMatrix should be the result returned by toMatrix().")
	
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
		if(dim(metagroupGenesMatrix)[2] != length(legendMg)) stop("The number of metagroups in the matrix and the legend do not match.")
	}
	
	# Initialize
	nMetagroups <- dim(metagroupGenesMatrix)[2]
	######## Set colors
	if(keepColors==TRUE)
	{
		colores <- setColors(as.character(sort(as.numeric(c(colnames(metagroupGenesMatrix), filteredOut)))))[colnames(metagroupGenesMatrix)]
	}else
	{
		colores <- setColors(colnames(metagroupGenesMatrix))
	}
	colnames(metagroupGenesMatrix) <- paste(grPrefix, colnames(metagroupGenesMatrix), sep=" ")

	#####################################################################################################
	#################################### Create Matrices  ###############################################	
	genesInManyMg <- which(apply(metagroupGenesMatrix, 1, sum)>1)
	
	mgJoinedGraph <- NULL
	if((length(genesInManyMg)>1) && any(plotType %in% c("static", "dynamic")))
	{
		mgJoined <- metagroupGenesMatrix[genesInManyMg,]
		
		# Replace 1 by mg number (for vertex color)
		for(col in 1:ncol(mgJoined))
		{
			thisCol <- which(mgJoined[,col]==1)
			mgJoined[thisCol,col] <- col
		}	
		
		if(!is.null(legendMg))
		{
			colnames(mgJoined) <- paste(colnames(mgJoined), legendMg, sep="\n")
		}
		
		# Create adjacency matrix (add nodes and mg to both, columns and rows)
		mgXmg <- matrix(data=0, ncol=dim(mgJoined)[2], nrow=dim(mgJoined)[2])
		rownames(mgXmg) <- colnames(mgJoined)
																
		genXgen <- matrix(data=0, ncol=dim(mgJoined)[1], nrow=dim(mgJoined)[1]+dim(mgJoined)[2])
		colnames(genXgen) <- rownames(mgJoined)
		
		mgJoined <- rbind(mgJoined, mgXmg)
		mgJoined <- cbind(genXgen, mgJoined)
		
		# Vertex colors (genes: white, mg:their color)
		vCols2 <- c(rep("white", length(genesInManyMg)), colores)
		
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
				tkplot(mgJoinedGraph, layout=vertexLayout, vertex.color=vCols2, edge.color=colores[E(mgJoinedGraph)$weight], edge.width=2, vertex.size=vSize, vertex.label.cex=vLabelCex)
			}else
			{
				if(plotType!= "none")
				{
					plot(mgJoinedGraph, layout=vertexLayout, vertex.color=vCols2, edge.color=colores[E(mgJoinedGraph)$weight], edge.width=2, vertex.size=vSize, vertex.label.cex=vLabelCex)
					title(main=plotTitle, sub=paste("Layout: ", layoutName, sep=""))
				}
			}
		}
	}
	
	if(returnGraph) return(mgJoinedGraph)
}













