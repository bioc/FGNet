# graph2plot:
# graphCommonGTsets = functionalNetwork (default)
# mgJoinedGraph = intersectionNetwork (bipartite)

plotFGNet <- function(graph2plot, plotType, plotOutput, plotExpression, vertexLayout, eWidth, eColor, vColor, vExpr, vSize, vLabelCex, colores, markGroup=NULL, trCols=NULL, plotTitle, plotTitleSub, plotLegend, legendPrefix, legendText)
{
    #####################################################################################################
    #########################################  Draw PLOTS  ##############################################
    if(plotOutput=="dynamic")
    {
        # if(is.null(plotIntersectionNetwork)) plotIntersectionNetwork <- TRUE
        if(ecount(graph2plot)>0)
        {
            vBorder <- "darkgrey"
            if(plotExpression == "fill") vBorder <- sub("#FFFFFF", "#888888", vColor)
            if(plotExpression == "border")
            {
                vBorder <- vExpr
            }
            
            tkplot(graph2plot, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label=V(graph2plot)$name, vertex.label.cex=vLabelCex, vertex.label.family="sans", vertex.color=vColor, vertex.size=vSize, vertex.frame.color=vBorder, vertex.label.color="black")  # vertex.label no estaba especificado para intesection network
        } else
        {
            warning("The dynamic plot is not available for networks without edges, try the 'static' plot.", immediate.=TRUE) # Tkplot error: Error in mapply(function(from, to, id) .tkplot.create.edge: zero-length inputs cannot be mixed with those of non-zero length
        }
    }
    if(plotOutput=="static")
    {        
        if(plotType=="bipartite")
        {
            add <- FALSE
            if(plotExpression == "border")
            {
                plot(graph2plot, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.color=vExpr, vertex.frame.color=vExpr,  vertex.label="",  vertex.size=ifelse(vExpr=="#888888",0,vSize+3)) 
                add <- TRUE
            }
            
        }else # "default" (FunctionalNetwork con background)
        {            
            add <- TRUE
            # Background and border
            if(plotExpression == "fill")
            {
                plot(graph2plot, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label="", mark.groups=markGroup, mark.col=trCols, mark.border=colores, vertex.size=0)
            }else
            { # Border or default/false
                plot(graph2plot, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label="", mark.groups=markGroup, mark.col=trCols, mark.border=colores, vertex.color=vExpr, vertex.frame.color=vExpr, vertex.size=ifelse(vExpr=="#888888",0,vSize+3))
            }
            
            # Legend
            if(plotLegend) 
            {
                if(!any(is.na(suppressWarnings(as.numeric(names(colores)))))) colorOrder <- order(as.numeric(names(colores)))
                else colorOrder <- order(names(colores))
                legend(-1.2, -1.2, legend=paste(legendPrefix, names(colores)[colorOrder], legendText, sep=""), fill=colores[colorOrder], bty="n", xjust=0, yjust=0)
            }
        }
        
        if(add) eWidth <-0
        plot(graph2plot, layout=vertexLayout, edge.width=eWidth, edge.color=eColor, vertex.label=V(graph2plot)$name, vertex.label.cex=vLabelCex, vertex.label.family="sans", vertex.color=vColor, vertex.size=vSize, vertex.frame.color=vExpr, vertex.label.color="black", add=add)   # en intersect. no estba vertex.label ni vertex.frame.color=vExpr,
        title(main=plotTitle, sub=plotTitleSub)    
    }
}