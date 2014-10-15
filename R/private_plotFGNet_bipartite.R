# metagroupsMatrix
# plotOutput

plotFGNet_bipartite <- function(metagroupsMatrix, keepAllNodes, plotAllMg, vLayout, plotOutput, plotExpression, eWidth, eColor, vExpr,vSize, vLabelCex, vColor, colores, plotTitle, plotTitleSub, plotLegend, legendPrefix, legendText)
{
    #####################################################################################################
    #################################### Create Matrices  ###############################################    
    genesInManyMg <- setNames(1:nrow(metagroupsMatrix), rownames(metagroupsMatrix))
    if(!keepAllNodes) genesInManyMg <- which(apply(metagroupsMatrix, 1, sum)>1)
    
    mgJoinedGraph <- NULL
    mgJoined <- NULL
    if(length(genesInManyMg)<=1)
    {
        warning("There is nothing to plot. There is no intersection between metagroups/clusters.", immediate.=TRUE)
    }else
    {
        if(!plotAllMg) metagroupsMatrix <- metagroupsMatrix[names(genesInManyMg),apply(metagroupsMatrix[names(genesInManyMg),,drop=FALSE], 2, sum)>0, drop=FALSE]
        
        mgJoined <- metagroupsMatrix[names(genesInManyMg),, drop=FALSE]
                        
        # Replace 1 by mg number (for vertex color)     
        if(!any(is.na(suppressWarnings(as.numeric(colnames(mgJoined))))))
        {
            weighted <- TRUE
            for(col in 1:ncol(mgJoined))
            {
                thisCol <- which(mgJoined[,col]==1)
                mgJoined[thisCol,col] <- colnames(mgJoined)[col]
            }      
        }else
        {
            weighted <- NULL
        }
        
        if(!is.null(legendText))
        {
            colnames(mgJoined) <- paste(colnames(mgJoined), legendText, sep="\n")
        }
        
        # Create adjacency matrix (add nodes and mg to both, columns and rows)
        mgXmg <- matrix(data=0, ncol=dim(mgJoined)[2], nrow=dim(mgJoined)[2])
        rownames(mgXmg) <- colnames(mgJoined)
        
        genXgen <- matrix(data=0, ncol=dim(mgJoined)[1], nrow=dim(mgJoined)[1]+dim(mgJoined)[2])
        colnames(genXgen) <- rownames(mgJoined)
        
        mgJoined <- rbind(mgJoined, mgXmg)
        mgJoined <- cbind(genXgen, mgJoined)
        
        # Vertex colors (genes: white, mg:their color)
        vColor <- setNames(c(rep("white", length(genesInManyMg)), colores[colnames(metagroupsMatrix)]), rownames(mgJoined))
        vShape <- ifelse(vColor=="white", "circle", "square")
        if(is.null(eWidth)) eWidth <- 2

        #####################################################################################################
        ####################################   Create Graph   ###############################################
        # Create graph
        mgJoinedGraph <- graph.adjacency(mgJoined, weighted=weighted, mode="undirected", diag=FALSE)

        #####################################################################################################
        #########################################  Draw PLOTS  ##############################################
        
        if(any(plotOutput %in% c("static", "dynamic")))
        {
            ###########################
            # Vertex layout
            if(is.null(vLayout)) vLayout<-"kk"
            if(!is.matrix(vLayout)) 
            {
                if(vLayout %in% "kk")
                {
                    vertexLayout <- layout.kamada.kawai(mgJoinedGraph)
                    layoutName <- "Kamada Kawai"
                }
                if(vLayout %in% "circle")
                {
                    vertexLayout <- layout.circle(mgJoinedGraph)
                    layoutName <- "Circle"
                }
                if(vLayout %in% "sugiyama") 
                { 
                    vertexLayout <- layout.sugiyama(mgJoinedGraph)$layout
                    layoutName <- "Sugiyama"
                }
            }else # Else: Layout is already provided
            {
                vertexLayout <- vLayout
                layoutName <- "Personalized"
            }
                  
            #######################
            # Edge colors
            eColor <- colores[as.character(E(mgJoinedGraph)$weight)]
            if(is.null(weighted))
            {
                eColor <- rep(NA, length(E(mgJoinedGraph)))
                for(vertName in colnames(metagroupsMatrix))
                {
                    vertInd <- which(V(mgJoinedGraph)$name == vertName)
                    eColor[incident(mgJoinedGraph, vertInd)] <- colores[vertName]
                }
            }
                        
            vExpr <- setNames(vExpr[names(vColor)], names(vColor))
            vExpr[which(is.na(vExpr))] <- "#888888"
            # Assign color only to "intersecton" nodes
            if(plotExpression=="fill") vColor[vColor=="white"] <- vExpr[vColor=="white"]
                        
            # Add legendPrefix
            if(legendPrefix!="")
            {
                clustIndex <- which(V(mgJoinedGraph)$name %in% colnames(metagroupsMatrix))
                V(mgJoinedGraph)$name[clustIndex] <- paste(legendPrefix, V(mgJoinedGraph)$name[clustIndex], sep="") # NO subir
            }
            #######################
            # Plot
            if(plotTitleSub=="") paste("Layout: ", layoutName, sep="")
            plotFGNet(graph2plot=mgJoinedGraph, plotType="bipartite", plotOutput=plotOutput, plotExpression=plotExpression, vertexLayout=vertexLayout, eWidth=eWidth, eColor=eColor, vColor=vColor,
                  vExpr=vExpr, vSize=vSize, vLabelCex=vLabelCex, vShape=vShape, colores=colores, markGroup=NULL, trCols=NULL,
                  plotTitle=plotTitle, plotTitleSub=plotTitleSub, plotLegend=plotLegend, legendPrefix=legendPrefix, legendText=legendText)
        }

        ### Return
        ret <- list(iGraph=list(bipartite=mgJoinedGraph), adjMat=list(bipartite=mgJoined))
        return(ret)
                    
    }# length(genesInManyMg)>1
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    