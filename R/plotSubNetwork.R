# txtFile <- "DavidClustering.txt"
# results <- getResults_david(txtFile)
# selCluster <- 1


### PENDIENTE:
# Gtlinker
# Top-GO... etc

plotSubNetwork <- function(results, selCluster, plotType=c("nw1","nw2","nw3","nw5","nw6","nw7","nw8"), plotOutput="dynamic", returnGraph=FALSE, geneExpr=NULL,plotExpression="border")
{
#   ret <- NULL
#   if(class(results)=="list")
#   {
#       grType <- names(results)[names(results)%in% c("clusters", "metagroups")]
#       if(grType == "clusters")
#       {
#         geneTermSets <- results$geneTermSets[results$geneTermSets$Cluster == selCluster,, drop=FALSE]
#         geneTermSets <- cbind(geneTermSets[,-1], Cluster=1:nrow(geneTermSets))
#         grType <- "Cluster"
#       }
#       if(grType == "metagroups")
#       {
#         geneTermSets <- results$geneTermSets[results$geneTermSets$Metagroup == selCluster,, drop=FALSE]
#         geneTermSets <- cbind(geneTermSets[,-1], Cluster=1:nrow(geneTermSets))
#         grType <- "Metagroup"
#       }
#   }else{
#       
#       
#       #### TEST
#       results <- cbind(Cluster=1:nrow(results),results)
#       
#       
#       
#       geneTermSets <- results[results$Cluster == selCluster,, drop=FALSE]
#       geneTermSets <- cbind(geneTermSets[,-1], Cluster=1:nrow(geneTermSets))
#       grType <- "Gene-term set"
#   }
#     
#   if(any(plotType %in% c("nw1","nw2","nw3")))
#   {
#     incidMat <- toMatrix(geneTermSets, key="genes")
#     if("nw1" %in% plotType) ret <- c(ret, functionalNetwork(incidMat, legendPrefix="gtSet ", plotTitle=paste(grType, selCluster,sep=""), plotTitleSub="Genes in gene-term sets", plotType=plotOutput, plotIntersectionNetwork=FALSE, returnGraph=returnGraph, returnAdjMatrices=returnGraph, geneExpr=geneExpr,plotExpression=plotExpression))  # Genes in GTsets
#     if("nw2" %in% plotType) ret <- c(ret, intersectionNetwork(incidMat, grPrefix="gtSet ", plotTitle=paste(grType, selCluster,sep=" "), plotTitleSub="Genes in gene-term sets", keepAllNodes = TRUE,plotType=plotOutput, returnGraph=returnGraph, geneExpr=geneExpr,plotExpression=plotExpression))
#     
#     if("nw3" %in% plotType) 
#     {
#       # DAVID - No intersection, GTLinker - YES.
#       incidMat <- toMatrix(geneTermSets, key="Terms")
#       ret <- c(ret, functionalNetwork(incidMat, legendPrefix="gtSet", plotTitle=paste(grType, selCluster,sep=" "),plotTitleSub="Terms in gene-term sets", plotIntersectionNetwork=FALSE, returnGraph=returnGraph, returnAdjMatrices=returnGraph, geneExpr=NULL)) # No intersection...  
#       # nw4? intersectionNetwork(incidMat, plotType="static", keepAllNodes=TRUE) (in david it does not work, there is only one term per gtset)
#     }
#   }
#   
#   
#   
#   ##### TO DO:  No vale para GTLINKER. Pensar
#   if(any(plotType %in% c("nw5","nw6","nw7","nw8")))
#   {    
#     rownames(geneTermSets) <- sapply(strsplit(geneTermSets$Terms, ":"), function(x) paste(capitalize(x[length(x)])," (" ,x[1], ")", sep=""))
#     # require(reshape2)
#     geneTermSetsReshaped <- melt(apply(geneTermSets[,c("Terms", "Genes")],1, function(x) unlist(strsplit(x["Genes"],","))))
#     geneTermSetsReshaped[,2] <- capitalize(sapply(strsplit(geneTermSetsReshaped[,2], ":"), function(x) x[length(x)]))
#     
#     colnames(geneTermSetsReshaped) <- c("Cluster","Terms")
#     matrices <- toMatrix(geneTermSetsReshaped, key="Terms")
#     # functionalNetwork(matrices, legendText=F, eColor="black") # Sin enlaces
#     if("nw5" %in% plotType) ret <- c(ret, functionalNetwork(matrices$clustersMatrix, legendText=FALSE ,plotTitle=paste(grType, selCluster,sep=" "), plotTitleSub="Terms sharing genes (?)", plotType=plotOutput, plotIntersectionNetwork=FALSE, returnGraph=returnGraph, returnAdjMatrices=returnGraph, geneExpr=NULL))   # TODOS con todos...
#     if("nw6" %in% plotType) ret <- c(ret, intersectionNetwork(matrices$clustersMatrix, plotTitle=paste(grType, selCluster,sep=" "), plotTitleSub="Terms - Genes", keepAllNodes=TRUE, plotType=plotOutput, returnGraph=returnGraph, geneExpr=NULL))
#     
#     
#     colnames(geneTermSetsReshaped)<- c("Genes","Cluster")
#     matrices <- toMatrix(geneTermSetsReshaped, key="Genes")
#     # functionalNetwork(matrices, legendText=F, eColor="black") # Sin enlaces
#     if("nw7" %in% plotType) ret <- c(ret, functionalNetwork(matrices$clustersMatrix, legendText=FALSE, plotTitle=paste(grType, selCluster,sep=" "), plotTitleSub="Genes sharing terms (?)", plotType=plotOutput, plotIntersectionNetwork=FALSE, returnGraph=returnGraph, returnAdjMatrices=returnGraph, geneExpr=NULL,plotExpression=plotExpression))   # TODOS con todos...
#     if("nw8" %in% plotType) ret <- c(ret, intersectionNetwork(matrices$clustersMatrix, plotTitle=paste(grType, selCluster,sep=" "), plotTitleSub="Genes - Terms", keepAllNodes=TRUE, plotType=plotOutput, returnGraph=returnGraph, geneExpr=geneExpr,plotExpression=plotExpression))
#   }
#   
#   
#   # BORRADOR DE GTLINKER:
# #   rownames(geneTermSets) <- geneTermSets$Terms
# #   mg12 <- melt(apply(geneTermSets[,c("Terms", "Genes")],1, function(x) unlist(strsplit(x["Genes"],","))))
# #   colnames(mg12)<- c("Metagroup","Terms")
# #   matrices <- toMatrix(mg12, key="Terms", clusterColumn="Metagroup")
# #   functionalNetwork(matrices$clustersMatrix, legendText=F, eColor="black", plotTitle=paste("Terms with common genes. Cl  ", cl,sep=""))#, plotType="dynamic")   # TODOS con todos...
# #   intersectionNetwork(matrices$clustersMatrix, plotTitle=paste("Terms - Genes. Cl  ", cl,sep=""), keepAllNodes=TRUE)#, plotType="dynamic")
# #   functionalNetwork(matrices, legendText=F, eColor="black") 
# #   
#   
#   
#   return(ret)
}  





# pdf("plots.pdf")
# plotSubNetwork(results, selCluster, plotOutput="static")
# dev.off()
