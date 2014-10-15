plotGoAncestors <- function(goIds, tColor=NULL, ontology=NULL, plotOutput="static", nCharTerm=50, nSize=NULL,labelCex=NULL, asp=NULL, fileName=NULL, height=1000)
{
    #### DB
    if(!is.list(ontology))
    {
        if(!loadInstPkg("GO.db")) stop("Package GO.db is required to plot GO trees.")
        
        # ontology == NULL: only OK if all GO ids are from the same DB 
        # ontology == character: selected DB is used
        # ontology == list: for internal use
        
        selectedDB <- ontology 
        if(!is.null(selectedDB))
        {
            if(!is.character(selectedDB) || length(selectedDB)>1) stop("ontology should be either 'BP', 'MF' *or* 'CC'")
            selectedDB <- toupper(selectedDB)
            if(!selectedDB %in% c("BP","MF","CC")) stop("ontology should be 'BP', 'MF' *or* 'CC'")
        }
        
        ontology <- list()
        ontology[["BP"]] <- AnnotationDbi::as.list(GOBPPARENTS)
        ontology[["MF"]] <- AnnotationDbi::as.list(GOMFPARENTS)
        ontology[["CC"]] <- AnnotationDbi::as.list(GOCCPARENTS)
        
        if(!is.null(selectedDB)) goIds <- goIds[goIds %in% names(ontology[[selectedDB]])]    
        if(is.null(selectedDB)) selectedDB <- names(which(sapply(ontology, function(x) any(goIds %in% names(x)))))
        if(length(selectedDB)!=1) stop("Provide 'ontology' or all GO IDs from the same type: BP, MF *or* CC.")
        ontology <- ontology[[selectedDB]]
    }
    
    adjList <- matrix(nrow=0, ncol=3)
    colnames(adjList) <- c("parent", "rel", "child")
    goChilds <- goIds
    
    # Find parents (calculate adjacency)
    while(any(goChilds != "all"))
    {
        goChilds <- goChilds[which(goChilds!="all")]
        
        goFamilyList <- lapply(ontology[goChilds], function(x) data.frame(parent=x, rel=names(x)))
        adjList <- rbind(adjList, as.matrix(cbind(do.call(rbind, goFamilyList), child=rep(names(goFamilyList),sapply(goFamilyList, nrow)))))
        
        # Next children
        goChilds <- unique(as.character(adjList[which(!adjList[,"parent"] %in% adjList[,"child"]),"parent"]))    
    }
    adjList <- adjList[which(adjList[,1]!="all"),, drop=FALSE]
    adjList <- adjList[,c("child","parent","rel"), drop=FALSE]
        
    #### Plot
    #adjList <- as.matrix(do.call(rbind.data.frame, adjList)) # unir
    goGraph <- graph.edgelist(adjList[,c("child","parent"), drop=FALSE], directed=TRUE)
    missingNodes <- rownames(goIds)[!goIds %in% unique(as.vector(adjList))]
    goGraph <- goGraph + vertex(missingNodes)
    if(vcount(goGraph)>0)
    {
        goGraph <- set.edge.attribute(goGraph, "type", value=adjList[,"rel"])
        
        BP <- "GO:0008150" %in% V(goGraph)$name
        goTerms <- sapply(V(goGraph)$name, function(x)
        {
            term <- capitalize(GOTERM[[x]]@Term)
            
            if (BP) # Shorten terms...
            {
                term <- gsub("( process)$", " pr.", term)
                term <- gsub(" development ", " dev. ", term)
                term <- gsub("( development)$", " dev.", term, perl=TRUE)
                term <- gsub("Development ", "Dev. ", term)
            }
            # Trim terms
            if(nchar(term)>nCharTerm) term <- paste(substr(term, 1, nCharTerm),"...", sep="")
            
            # split in two lines
            #gsub("(.*) (.*) (.*)", "\\1 \n\\2 \\3", term)
            spaceLoc <- gregexpr(" ", term, fixed=TRUE)[[1]]
            if((spaceLoc[1]!=-1) && (length(spaceLoc)>1))
            {
                spaceLoc <- spaceLoc[ceiling(length(spaceLoc)/2)]
                term <- paste(substr(term, 1, spaceLoc-1), substr(term, spaceLoc+1, nchar(term)), sep="\n")
            }

            term
        })
        names(goTerms) <- NULL
        goGraph <- set.vertex.attribute(goGraph, "term", value=goTerms) # comprobar
        
        # Node color
        goGraph <- set.vertex.attribute(goGraph, "initial", value=V(goGraph)$name%in% goIds)
        finalLeaves <- igraph::degree(goGraph,mode="in")==0
        
        if(plotOutput!="none")
        {
            vColor <- c("slategray1", "skyblue2") # Default: color for initial terms, darker if leaf.
            if(length(tColor) == 1) vColor <- rep(tColor,2) # If only one color is provided, all initial terms are coloured with it        
            # Apply to graph...
            goGraph <- set.vertex.attribute(goGraph, "color", value=c("white",vColor[1])[as.numeric(V(goGraph)$initial)+1]) # comprobar
            V(goGraph)$color[finalLeaves] <- vColor[2]
            
            if(!is.null(names(tColor))) V(goGraph)$color <- tColor[V(goGraph)$name]
        
            # Calculate node level for layout
            tmpGoGraph <- goGraph
            nodeLevels <- NULL
            while(any((leaves <- igraph::degree(tmpGoGraph,mode="in"))==0))
            {
                leaves <- leaves[leaves==0]
                nodeLevels <- c(nodeLevels, leaves)
                nodeLevels <- nodeLevels +1
                
                tmpGoGraph <- delete.vertices(tmpGoGraph, names(leaves))
            }
            
            # Layout
            #nodeLevels <- nodeLevels + sample(c(-0.1,+0.1), length(nodeLevels), replace=TRUE)
            vertexLayout <- layout.sugiyama(goGraph, layers=nodeLevels[V(goGraph)$name], hgap=5,maxiter=300, attributes="all")
            eColors <- c("darkgray","cadetblue3", "darkgoldenrod1","pink", "darkolivegreen1")
            names(eColors) <- c("is_a","part_of","regulates","positively_regulates","negatively_regulates")
            dColorsSimplif <- c("gray", "blue", "yellow","pink", "green") # Simplified name        
            if(length(unique(E(goGraph)$type[!E(goGraph)$type %in% names(eColors)]))>0) warning("missing rel")
            
            if(is.null(labelCex)) labelCex <-  ifelse(vcount(goGraph)<20, 1, 0.75)
        
            # if... plotDynamic
            if(plotOutput=="dynamic")
            {
                if(is.null(labelCex)) labelCex <- 1
                if((apply(vertexLayout$layout,2, sum)[1]==0) && TRUE)  vertexLayout$layout[1,1] <- vertexLayout$layout[1,1]+0.1
                tkplot(goGraph, layout=vertexLayout$layout, 
                             edge.color=eColors[E(goGraph)$type], edge.arrow.size=0.8,
                             vertex.color=V(goGraph)$color, vertex.frame.color=V(goGraph)$color,
                             vertex.label=V(goGraph)$term, vertex.label.color="black",  vertex.label.cex=labelCex,
                             vertex.size=15)
                print(paste("Links: ", paste(names(eColors), " (", dColorsSimplif,")", collapse=", ", sep=""), collapse=""))
            } else # Not dynamic
            {        
                if(!is.null(fileName))
                {
                    if(is.null(asp))
                    {
                        asp <- ((max(vertexLayout$layout[,2])-min(vertexLayout$layout[,2]))/max(abs(vertexLayout$layout[,1])))*3
                        asp <- min(asp, 1)
                    }
                    if(is.null(nSize)) nSize <- 15*asp
                    if(is.null(labelCex)) labelCex <-  1+(1*asp)                    
                    
                    png(fileName, height = height, width = height/asp)
                }
                if(is.null(nSize)) nSize <- ifelse(vcount(goGraph)<20, 15, 10)
                
                arrowSize <- nSize/30
                
                origvert <- c(rep(TRUE, vcount(goGraph)), rep(FALSE, nrow(rbind(vertexLayout$layout.dummy))))
                realedge <- !is.na(get.edgelist(vertexLayout$extd_graph)[,2])        
                
                
                plot(vertexLayout$extd_graph, #layout=vertexLayout$layout,
                         edge.color=eColors[E(vertexLayout$extd_graph)$type],
                         edge.arrow.mode=ifelse(realedge, 2, 0), edge.arrow.size=arrowSize,
                         vertex.shape=ifelse(origvert, "circle", "none"),
                         vertex.size=ifelse(origvert, nSize, 0),
                         vertex.color=V(vertexLayout$extd_graph)$color, vertex.frame.color=V(vertexLayout$extd_graph)$color,
                         vertex.label=ifelse(origvert, V(goGraph)$term, ""), vertex.label.color="black",  vertex.label.cex=labelCex,
                         asp=asp,margin=0)        
                
                legend(-1.2,1.2,title="Link color:",legend = levels(factor(E(goGraph)$type)), col = eColors[levels(factor(E(goGraph)$type))], lty = 1, lwd=3, bg="transparent") 
                if(!is.null(fileName)) dev.off()
            }
        }
        
        finalLeaves <- cbind(sort(sapply(V(goGraph)$name[finalLeaves], function(x) capitalize(GOTERM[[x]]@Term))))
        finalLeaves <- cbind(rownames(finalLeaves), finalLeaves)
    }
   
    
    ret <- list(iGraph=goGraph, leaves=finalLeaves)
    invisible(ret)
}

