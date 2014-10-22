
# incidMatrices <- fea2incidMat(res_topGo)
analyzeNetwork <- function(incidMatrices, fNw=NULL, plotOutput=TRUE, colors=NULL)
{
    if(is.null(fNw)) fNw <-functionalNetwork(incidMatrices, plotOutput="none")
     
    clNw <- graph.adjacency(fNw$adjMat$commonClusters, weighted=TRUE, mode="undirected", diag=FALSE)  # weigted: only one edge between genes (returned is NOT weighted)
    gtsNw <- fNw$iGraph$commonGtSets
    
    data("groupTypes", envir = environment())
    groupTypes<- get("groupTypes", envir  = environment())
    clMatrix <- incidMatrices[[which(names(incidMatrices)%in%groupTypes[,"tablesName"])]]
    
    #########################################
    # Global betweenness:
    glBetw <- round(sort(betweenness(clNw), decreasing=TRUE))
    # Global degree
    glDegr <- round(sort((igraph::degree(clNw)/(length(V(clNw))-1)*100), decreasing=TRUE),digits=2)
    # Equivalent to boxplot: cbind(degree=round(summary(glDegr))[-4], betweenness=round(summary(glBetw))[-4])
    
    
    ## Cluster betweenness (intra-modular hub)
    genesMg <- apply(clMatrix,2, function(x) names(x)[which(x==1)])
    clSubgraphs <- list()
    clBetw <- list(Global=glBetw)
    clDegr <- list(Global=glDegr)
    for(cl in names(genesMg))
    {
        clSubgraphs[[cl]] <- induced.subgraph(gtsNw, which(V(gtsNw)$name %in% genesMg[[cl]]))
        clBetw[[cl]] <- round(sort(betweenness(clSubgraphs[[cl]]), decreasing=TRUE))
        clDegr[[cl]] <- round(sort(igraph::degree(clSubgraphs[[cl]]), decreasing=TRUE)/(length(V(clSubgraphs[[cl]]))-1),digits=2)*100
    }    
    
    #########################################
    ## Plot
    if(plotOutput)
    {
        if(is.null(colors)) colors <- rep("white", length(clDegr))
        par(mfrow=c(1,2))
        boxplot(clDegr, main="Normalized node degree", sub="Global / intra-cluster", outpch=16,outcol=c("black", colors), outlwd=1, boxfill=c("black", colors), medcol=c("black", colors), medlwd = 10, lwd = 1, ylim=c(0,100), ylab="Percentage of total nodes", axes=FALSE)
        box()
        axis(side=1, labels=names(clDegr), at=1:length(clDegr), las=2)
        axis(side=2, las=1)
        
        boxplot(c(0,clBetw[-1]), main="Node betweenness", sub="Global / intra-cluster", outpch=16, outcol=c("white", colors),outlwd=1, boxfill=c("white", colors), medcol=c("white", colors), medlwd = 10, lwd = 1, axes=FALSE)
        maxBetwGlobal <- max(clBetw$Global)
        if(maxBetwGlobal>0) 
        {
            transf <- max(unlist(clBetw[-1])/maxBetwGlobal)
            betwGlobal <- clBetw$Global*transf
        }else
        {
            betwGlobal <- 0
        }
        boxplot(betwGlobal, add=TRUE, outpch=16, outcol="black", outlwd=1, boxfill="black", medcol="black", medlwd = 10, lwd = 1, axes=FALSE)
        box()
        abline(v=1.5)
        axis(side=1, labels=names(clDegr), at=1:length(clDegr), las=2)
        axis(side=4, las=1)
        axis(side=2, las=1, labels=c(0, max(clBetw$Global)), at=c(0, max(unlist(clBetw[-1])))) # Global side
    }
    
    #########################################
    ## Stats
    ret <- list()
    
    # Node degree and betweenness:
    ret[["degree"]] <- clDegr
    ret[["betweenness"]] <- clBetw
    
    # Clustering coeficient. Transitivity measures the probability that the adjacent vertices of a vertex are connected.
    ret[["transitivity"]] <- c(commonClustersNw=transitivity(clNw, type="undirected"), commonGtSetsNw=transitivity(gtsNw, type="undirected"))
    
    # Hubs
    clBetw <- sapply(clBetw, function(x) x[x!=0]) # Remove genes with betw==0
    clBetw <- clBetw[sapply(clBetw[names(clBetw)!="Global"], length)>0] # Remove clusters without genes with betw>0 (except Global)
    
    hubs <- unique(unlist(sapply(clBetw, names))) # Gene names
    hubs <- matrix(0, nrow=length(hubs), ncol=length(clBetw), dimnames=list(hubs, names(clBetw)))
    for(cl in names(clBetw))
    {
        hubs[,cl] <- clBetw[[cl]][rownames(hubs)]
    }
    hubs[which(is.na(hubs))] <- 0
    ret[["betweenessMatrix"]] <- hubs
    
    hubsList <- lapply(clBetw, function(x){
       names(x[x>=quantile(x,probs=0.75)])
    })
    ret[["hubsList"]] <- hubsList
    
    # Intra-modular hubs count
    intraModularHubs <- hubsList[names(hubsList)!="Global"]
    tmpHubsTable <- table(unlist(intraModularHubs))
    tmpHubsTable <- sort(tmpHubsTable[tmpHubsTable>1], decreasing=TRUE)
    ret[["intraHubsCount"]] <- tmpHubsTable
    
    invisible(ret)    
}



# 
# #########################################
# ### Centralization:
# centralization.degree (clNw, mode ="total")
# # $centralization
# # [1] 0.599435
# # 
# # $theoretical_max
# # [1] 3540
# 
# centralization.closeness (clNw, mode ="total")
# # $centralization
# # [1] 0.4355989
# # 
# # $theoretical_max
# # [1] 29.24786
# 
# centralization.betweenness (clNw, directed = FALSE)
# # $centralization
# # [1] 0.08470773
# # 
# # $theoretical_max
# # [1] 100949
# 
# centralization.evcent (clNw, directed = FALSE)
# # $centralization
# # [1] 0.6417392
# # 
# # $theoretical_max
# # [1] 58
