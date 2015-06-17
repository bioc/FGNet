
# Version 2.1, ADDED: geneExpr=NULL,plotExpression=c("border","fill")
# Quitar: metagroupAttribute=NULL --> columnas en metagroupsMatrix: ordenadas por el atributo
#   vertexLayout:  previous Vlayout (para todos) o  "kk", "circle"...etc: Solo se usa en intersection. 


# functionalNetwork:
# bgTransparency=0.4
# eColor="#323232"
# weighted=FALSE
# 
# intersectionNetwork 
# keepAllNodes=FALSE     # if(!keepAllNodes) genesInManyMg <- which(apply(metagroupsMatrix, 1, sum)>1)
# plotAllMg=FALSE


functionalNetwork <- function(incidMatrices, plotType=c("default", "bipartite")[1], plotOutput="static", plotTitle="Functional Network", plotTitleSub=NULL, legendPrefix=NULL, legendText=NULL,  geneExpr=NULL, plotExpression=c("border","fill"), vExprColors=c(neg="#008000", zero="white", pos="#FF2020"), vSize=12, vLabelCex=2/3, vLayout=NULL,keepColors=TRUE, bgTransparency=0.4, eColor="#323232", eWidth=NULL, weighted=FALSE, keepAllNodes=FALSE, plotAllMg=FALSE)
{    
    #####################################################################################################
    ########################### Check & initialize arguments ############################################
    filteredOut <- NULL
    gtSetsMatrix <- NULL
    
    if(!any(plotType %in% c("default", "bipartite"))) stop('plotType should be "default" or "bipartite.')
    # incidMatrices
     grType <- "Cluster"
     if(is.list(incidMatrices) && ("metagroupsMatrix" %in% names(incidMatrices)))   grType <- "Metagroup"
    if(is.matrix(incidMatrices)) metagroupsMatrix <- incidMatrices
    if(is.list(incidMatrices) && any(c("metagroupsMatrix", "clustersMatrix") %in% names(incidMatrices)))
    {
        if("metagroupsMatrix" %in% names(incidMatrices)) metagroupsMatrix <- incidMatrices$metagroupsMatrix
        if("clustersMatrix" %in% names(incidMatrices)) metagroupsMatrix <- incidMatrices$clustersMatrix
        if("gtSetsMatrix" %in% names(incidMatrices))  gtSetsMatrix <- incidMatrices$gtSetsMatrix
        if("filteredOut" %in% names(incidMatrices))   filteredOut <- incidMatrices$filteredOut
    }
    if(is.matrix(metagroupsMatrix) && is.null(gtSetsMatrix))
    {
        gtSetsMatrix <- metagroupsMatrix
        filteredOut <- NULL
        
        #warning("Only metagroupsMatrix provided.")
        if(plotTitle=="Functional Network") plotTitle <- "Metagroups Network"
    }    
    if(!is.matrix(metagroupsMatrix))  stop("metagroupsMatrix not valid. 'incidMatrices' should be the result returned by toMatrix().")
    if(!is.matrix(gtSetsMatrix))  stop("gtSetsMatrix not valid. 'incidMatrices' should be the result returned by toMatrix().")

    # plotOutput
    if(!is.character(plotOutput))  stop('plotOutput should be either "static", "dynamic" or "none".') 
    plotOutput <- tolower(plotOutput)
    if(!plotOutput %in% c("static", "dynamic", "none")) stop('plotOutput should be in "static", "dynamic" or "none".') 

    # plotTitle & plotTitleSub
    if(!is.character(plotTitle)) stop("plotTitle should be character.")
    if(!is.null(plotTitleSub) && !is.character(plotTitleSub)) stop("plotTitle should be character.")
    if(is.null(plotTitleSub))
    {
        plotTitleSub <- setNames(c("","Intersection network"), c("default","bipartite"))[plotType]     ### TO DO? (if keep all nodes... also intersection??)
    }
    
    # legendText & legendPrefix
    plotLegend <- TRUE
    if(is.logical(legendText) && legendText==FALSE)
    {
        plotLegend <- FALSE
        legendText <- NULL
    }
    if(!is.null(legendText))
    {
        if(!is.character(legendText)) stop("legendText should be a character vector containing the description of each cluster.")
        if(ncol(metagroupsMatrix) != length(legendText)) stop("The number of clusters/metagroups in the matrix and the legend do not match.")
        if(plotType=="default") legendText <- paste(": ", legendText, sep="")
    }    
    
    if(is.null(legendPrefix))
    {
        data("groupTypes", envir = environment())
        groupTypes<- get("groupTypes", envir  = environment())
        legendPrefix <- groupTypes[grType, "prefix"]
    }
    if(!is.character(legendPrefix)) stop("legendPrefix should be character.")
    
    ## Expression
    if(!is.null(geneExpr) && !is.numeric(geneExpr)) stop('geneExpr should be a numeric vector.')
    if(!is.character(plotExpression)) stop("plotExpression should be either 'border' or 'fill'.")
    plotExpression <- tolower(plotExpression[1])
    if(!plotExpression %in% c("border","fill")) stop("plotExpression should be either 'border' or 'fill'.")
    if(length(vExprColors)!=3) stop("vExprColors should be a character vector of length 3: colors for negative, zero and positive values.")
    if(is.null(names(vExprColors))) names(vExprColors) <- c("neg", "zero", "pos")
    # Other plot options
    if(!is.numeric(vSize))  stop("vSize should be numeric.")
    if(!is.numeric(vLabelCex))  stop("vLabelCex should be numeric.")
    if(!is.logical(keepColors)) stop("keepColors should be TRUE or FALSE")
    
    # Specific options
    if("default" %in% plotType)
    {
        if(!is.null(bgTransparency) && !is.numeric(bgTransparency))  stop("bgTransparency should be numeric.")
        if(!is.logical(weighted)) stop("weighted should be either TRUE or FALSE.")
        
        if(!is.null(vLayout) && (!is.matrix(vLayout) || (dim(vLayout)[2] != 2))) stop("Non-valid layout.")
    }
    
    #####################################################################################################
    ####################################  Plot settings  ################################################
    nGenes <- dim(metagroupsMatrix)[1]
    nMetagroups <- dim(metagroupsMatrix)[2]
    
    # In case of term nodes keep only their description as label
    if(grepl(":", rownames(metagroupsMatrix)[1]))
    {
        clMgMxNames <- clDescriptions(rownames(metagroupsMatrix))
        rownames(metagroupsMatrix) <- paste(clMgMxNames[,"TermDescription"], " (", clMgMxNames[,"Annotation"], ")", sep="")
        
        clGtSetMxNames <- clDescriptions(rownames(gtSetsMatrix))
        rownames(gtSetsMatrix) <- paste(clGtSetMxNames[,"TermDescription"], " (", clGtSetMxNames[,"Annotation"], ")", sep="")
    }
    
    ######## Set cluster color
    if(any(is.na(suppressWarnings(as.numeric(filteredOut))))) keepColors<- FALSE
    if(keepColors==TRUE)
    {
        if(all(!is.na(suppressWarnings(as.numeric(colnames(metagroupsMatrix))))))
        {
            allMgSorted <- as.character(sort(as.numeric(c(colnames(metagroupsMatrix), filteredOut))))
            allMgNotSorted <- as.character(as.numeric(colnames(metagroupsMatrix)))
        } else {
            allMgSorted <- sort(as.character(c(colnames(metagroupsMatrix), filteredOut)))
            allMgNotSorted <- as.character(colnames(metagroupsMatrix))
        }
        colores <- setColors(allMgSorted)[allMgNotSorted]
        trCols <- setColors(allMgSorted, transparency=bgTransparency)[allMgNotSorted]
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
    if(!is.list(markGroup))
    {
        if(!is.matrix(markGroup))  # in case of only one gene/term per metagroup
        {
            markGroup <- t(as.matrix(markGroup)) #  if(ncol(metagroupsMatrix)>1) ?
        }        
        markGroup <- split(markGroup, rep(1:ncol(markGroup), each = nrow(markGroup))) #in case of same number of vertex in each group...
    }
    
    #######################################################
    # Set color for expression
    if(!is.null(geneExpr))
    {
        geneExpr <- geneExpr[rownames(gtSetsMatrix)] # Subset, in case the whole eset was provided
        
        vExpr <- rep("#888888", length(geneExpr)) # Grey
        names(vExpr) <- names(geneExpr)  

        if(plotExpression == "fill")
        {
            vExprColorsRGB <- col2rgb(vExprColors, alpha = FALSE)/255
            vExpr[which(geneExpr == 0)] <- "white"
#             if(any(geneExpr < 0,na.rm=TRUE))   vExpr[which(geneExpr < 0)]  <- color.scale( geneExpr[which(geneExpr < 0)], cs1=c(0,1), cs2=c(0.5,1),cs3=c(0,1),alpha=1, xrange=range(min(geneExpr,na.rm=TRUE),0))
#             if(any(geneExpr > 0,na.rm=TRUE))   vExpr[which(geneExpr > 0)]  <- color.scale( geneExpr[which(geneExpr > 0)], cs1=c(1,1), cs2=c(1,0),cs3=c(1,0), alpha=1, xrange=range(0, max(geneExpr,na.rm=TRUE)))
            if(any(geneExpr < 0,na.rm=TRUE))   
            {
                if(length(table(geneExpr[which(geneExpr < 0)]))==1)
                {
                    vExpr[which(geneExpr < 0)] <- vExprColors["neg"]            
                }else
                {
                    vExpr[which(geneExpr < 0)] <- color.scale(geneExpr[which(geneExpr < 0)], cs1=vExprColorsRGB["red",c("neg","zero")], cs2=vExprColorsRGB["green",c("neg","zero")], cs3=vExprColorsRGB["blue",c("neg","zero")], alpha=1, xrange=range(min(geneExpr,na.rm=TRUE),0))                   
                }
            }
            if(any(geneExpr > 0,na.rm=TRUE))  
            {
                if(length(table(geneExpr[which(geneExpr > 0)]))==1)
                {
                    vExpr[which(geneExpr > 0)] <- vExprColors["pos"]            
                }else
                {
                    vExpr[which(geneExpr > 0)]  <- color.scale( geneExpr[which(geneExpr > 0)], cs1=vExprColorsRGB["red",c("zero","pos")], cs2=vExprColorsRGB["green",c("zero","pos")], cs3=vExprColorsRGB["blue",c("zero","pos")], alpha=1, xrange=range(0, max(geneExpr,na.rm=TRUE)))
                }
            }
                
            vColor <- vExpr
        }
        if(plotExpression == "border")
        {
            # Expresion en borde de nodo. 
            ##### Borde aparte (2 plots: back-expresion, front-mg)
            vExpr[which(geneExpr == 0)]  <- vExprColors["zero"]
            vExpr[which(geneExpr > 0)]  <-  vExprColors["pos"] # "#FF4949"
            vExpr[which(geneExpr < 0)]  <-  vExprColors["neg"] #"#008000"
        }
    } else { 
        vExpr <- "#888888" 
        plotExpression <- FALSE        
    }

print(vExpr)
    
## RECOLOCAR:     geneExpr[]
# PARA default: V(graphCommonGTsets)$name
# Para bipartite: V(....)$name    

    
    
    ########################################
    # Plots: 
    ret <- NULL
    if("default" %in% plotType)
    {
        ret <- plotFGNet_default(metagroupsMatrix, gtSetsMatrix, vLayout=vLayout, plotOutput=plotOutput, plotExpression=plotExpression, eWidth=eWidth, eColor=eColor, vColor=vColor,
                            vExpr=vExpr, vSize=vSize, vLabelCex=vLabelCex, colores=colores,markGroup=markGroup, trCols=trCols, weighted=weighted,
                            plotTitle=plotTitle, plotTitleSub=plotTitleSub, plotLegend=plotLegend, legendPrefix=legendPrefix, legendText=legendText)
    }
    
    if("bipartite" %in% plotType)
    {
        ret_bipartite <- plotFGNet_bipartite(metagroupsMatrix=metagroupsMatrix, keepAllNodes=keepAllNodes,plotAllMg=plotAllMg, vLayout=vLayout, plotOutput=plotOutput, plotExpression=plotExpression, eWidth=eWidth, eColor=eColor,  vColor=vColor,
                  vExpr=vExpr, vSize=vSize, vLabelCex=vLabelCex, colores=colores, 
                  plotTitle=plotTitle, plotTitleSub=plotTitleSub, plotLegend=plotLegend, legendPrefix=legendPrefix, legendText=legendText)
        
        if(is.null(ret)) # Merge with default
        {
            ret <- ret_bipartite
        }else
        {
            ret$iGraph <- c(ret$iGraph, ret_bipartite$iGraph)
            ret$adjMat <- c(ret$adjMat, ret_bipartite$adjMat)
        }
    }
    invisible(ret)
}


