# 
# FGNET - DEFAULT (genes-genes with background)
# 
plotFGNet_default <- function(metagroupsMatrix,gtSetsMatrix, vLayout, plotOutput, plotExpression, eWidth, eColor, vColor, vExpr, vSize, vLabelCex, colores, markGroup, trCols, weighted, plotTitle, plotTitleSub, plotLegend, legendPrefix, legendText)
{
    #####################################################################################################
    #################################### Create Matrices  ###############################################
    nGenes <- nrow(metagroupsMatrix)
    nMetagroups <- ncol(metagroupsMatrix)
    # Gene - Gene edges 
    nCommonGTsets <- matrix(ncol=nGenes, nrow=nGenes, data=0)
    rownames(nCommonGTsets) <- rownames(metagroupsMatrix)
    colnames(nCommonGTsets) <- rownames(metagroupsMatrix)
    # Counts
    nCommonMgroups <- nCommonGTsets
    # Fill:
    for(gen1 in rownames(metagroupsMatrix))
        for(gen2 in rownames(metagroupsMatrix))
        {
            if(gen1!=gen2) nCommonMgroups[gen1, gen2] <- sum(metagroupsMatrix[gen1,] + metagroupsMatrix[gen2,]==2)
            if(gen1!=gen2) nCommonGTsets[gen1, gen2] <- sum(gtSetsMatrix[gen1,] + gtSetsMatrix[gen2,]==2)
        }
    
    #####################################################################################################
    ####################################   Create Graph   ###############################################
    # weighted=TRUE The name of the edge attribute will be weight.
    # weighted=NULL (Unweighted graph. The elements of the adjacency matrix give the number of edges between vertices.)
    graphCommonGTsets <- graph.adjacency(nCommonGTsets, weighted=TRUE, mode="undirected", diag=FALSE) 
    graphCommonMgroups <- graph.adjacency(nCommonMgroups, weighted=NULL, mode="undirected", diag=FALSE)
    
    #####################################################################################################
    #########################################  Draw PLOTS  ##############################################
    if(any(plotOutput %in% c("static", "dynamic")))
    {
        ###########################
        # Vertex layout
        # weighted=NULL (Unweighted graph. The elements of the adjacency matrix give the number of edges between vertices.) -> No parece importar para el layout
        vertexLayout <- vLayout
        if(is.null(vertexLayout))
        {
            # Set layout based on the metagroups
            if(nMetagroups>1) 
            {
                # (graphCommonMgroups) NOT weigted: place genes with more mg in common closer
                graph4layout <- graphCommonMgroups
            } else { # If there is only one metagroup, set layout based on the common gene-term sets
                graphCommonMgroups <- NA
                graph4layout <- graphCommonGTsets   
            }        
            
            if(FALSE)
            {
                # Parecido a fruchterman reingold, pero apelotona un poco menos los nodos comunes
                vertexLayout <- layout.kamada.kawai(graph4layout)    
            }else
            {
                # Si son dos grupos separados (no tienen ningun nodo en comun) kamada kawai superpone los dos grupos = MALO
                vertexLayout <- layout.fruchterman.reingold(graph4layout)#,  area=200000, repulserad=200000*vcount(graph4layout))
            }
        }
        
        #######################
        # Edge weight
        if(is.null(eWidth))
        {
            if(weighted && ecount(graphCommonGTsets)>0)
            {
                eWidth <- E(graphCommonGTsets)$weight
                eWidth <- floor((eWidth/max(eWidth))*5)
                if(min(eWidth)<1) eWidth <- eWidth +1
                if(min(eWidth) == max(eWidth)) eWidth <- 1
            }else eWidth <- 1
        }
        
        #######################
        # Plot
        plotFGNet(graph2plot=graphCommonGTsets, plotType="default", plotOutput=plotOutput, plotExpression=plotExpression, vertexLayout=vertexLayout, eWidth=eWidth, eColor=eColor, vColor=vColor,
                  vExpr=vExpr, vSize=vSize, vLabelCex=vLabelCex, colores=colores, markGroup=markGroup, trCols=trCols,
                  plotTitle=plotTitle, plotTitleSub=plotTitleSub, plotLegend=plotLegend, legendPrefix=legendPrefix, legendText=legendText)   
    }

    ret <- list(iGraph=list(commonClusters=graphCommonMgroups, commonGtSets=graphCommonGTsets),
            adjMat=list(commonClusters=nCommonMgroups, commonGtSets=nCommonGTsets))

    return(ret)
}