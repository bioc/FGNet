# Functions:
#     buildTermsTable
#     goTreeLinks
#     createHtml - MAIN

# (Create html links)
buildTermsTable <- function(termsDescriptions, leaves, boldTerms)
{
    if(!is.null(leaves))
    {
        if(grepl("Cluster", names(termsDescriptions)[[1]]))    names(leaves) <- paste("Cluster",names(leaves))
        if(grepl("Metagroup", names(termsDescriptions)[[1]]))    names(leaves) <- paste("Metagroup",names(leaves))
        
        for(n in names(termsDescriptions))
        {
            gos <- grepl("GO", termsDescriptions[[n]][,"TermID"])
            
            if(any(gos)) termsDescriptions[[n]] <- rbind(termsDescriptions[[n]][!gos,],termsDescriptions[[n]][gos,, drop=FALSE][termsDescriptions[[n]][gos,"TermID", drop=FALSE] %in% leaves[[n]],])
        }
    }
    
    termLinks <- lapply(termsDescriptions, function(mgMx) {
        if(!is.matrix(mgMx)) mgMx <- t(as.matrix(mgMx)) # 1 row?
        
        if(!is.null(boldTerms))
        {
            bold <- as.character(boldTerms[boldTerms %in% rownames(mgMx)])
            mgMx[bold,"TermDescription"] <- paste("<b>",mgMx[bold,"TermDescription"],"</b>",sep="") 
        }
        
        descrLinks <- apply(mgMx, 1, function(term){
            
            if(is.na(term["Link"]))
            {
                link <- term["TermDescription"]
            }else
            {
                th1 <- th2 <- th3 <- ""
                if(grepl(".png", term["Link"]))  ## Add thumbnail?
                {
                    th1 <- "<div id='thumbnail'>"
                    th2 <- paste("<img src='",term["Link"],"'/>", sep="")
                    th3 <- "</div>"
                }   
                link <- paste(th1, "<a href='", term["Link"],"' target='_blank'>", term["TermDescription"] ,th2, "</a>", th3, sep="")                
            }
            return (link)
        })
        
        descrLinks <- cbind(descrLinks, mgMx[,"Annotation"])
        colnames(descrLinks) <- NULL
        rownames(descrLinks) <- NULL
        return(descrLinks)
    })
    
    return(list(termLinks=termLinks, noLinks=lapply(termsDescriptions, function(x) x[,"TermDescription",drop=FALSE])))
}

# Creates the GO tree and returns the link for each metagroup
goTreeLinks <- function(gtSets, grPrefix, geneExpr, folder, plotGoTree, clusterColumn)
{   
    if(!"GO.db" %in% rownames(installed.packages()))
        {
            tryCatch({
                    source("http://bioconductor.org/biocLite.R")
                    biocLite(GO.db)
            }, error = function(e) {})
    } 
    require(GO.db)
    
    data("groupTypes", envir = environment())
    groupTypes<- get("groupTypes", envir  = environment())
        
    # gtSets <- notFilteredGeneTermSets
    fileNames <- list()
    leaves <- list()
    
    # Initialize gtSets     
    if(grPrefix=="mg") 
    {
        gtSets  <- do.call(rbind,apply(as.matrix(gtSets),1, function(x){
            cbind(x["Metagroup"], strsplit(x["Terms"],";")$Terms, x["Genes"])
        }))
        colnames(gtSets) <- c("Metagroup","Term","Genes")
        rownames(gtSets) <- NULL
        gtSets <- data.frame(gtSets)
    }
    
    # Keep only GO terms
    gtSets <- gtSets[grepl(pattern="GO:",gtSets$Terms),]
    
    if(nrow(gtSets)>0 && plotGoTree)
    {
        # Create DB
        parentsDB <- list()
        parentsDB[["BP"]] <- AnnotationDbi::as.list(GOBPPARENTS)
        parentsDB[["MF"]] <- AnnotationDbi::as.list(GOMFPARENTS)
        parentsDB[["CC"]] <- AnnotationDbi::as.list(GOCCPARENTS)
        
        # For each Mg/Cluster...
        geneSplit <- ","
        clusterColumn <- rownames(groupTypes)[which(groupTypes[,"prefix"]==grPrefix)]
        for(cl in as.character(unique(gtSets[[clusterColumn]])))
        {
            gtsetsCl <- gtSets[gtSets[,clusterColumn]==cl,]
            # Gtlinker: Merge genes from same term in different gtset
            if(grPrefix=="mg")
            {
                tmp <- cbind(as.character(unique(gtsetsCl$Terms)), NA)
                colnames(tmp) <- c("Term", "Genes")
                tmp[,"Genes"] <- apply(tmp, 1, function(x) paste(unique(unlist(lapply(gtsetsCl[gtsetsCl$Terms==x[1],"Genes"], function(x) strsplit(as.character(x),split=geneSplit)[[1]]))), collapse=geneSplit))
                gtsetsCl <- data.frame(tmp)
            }     
            
            # Count terms UP and DOWN (para gtlinker se podria optimizar...)
            if(is.null(geneExpr))
            {
                termExprs <- NULL
            }else{ 
                termExprs <- sapply(gtsetsCl$Genes, function(x) {
                                    genes <- strsplit(as.character(x),split=geneSplit)[[1]]
                                    (sum(geneExpr[genes]>0)-sum(geneExpr[genes]<0))/length(genes) # NA if any missing...
                                    })
            }
            
            # Set term color according to geneExpr
            termsAll <- as.character(gtsetsCl$Terms)
            termsAll <- sapply(strsplit(termsAll,":"), function(x) paste(x[1],x[2],sep=":"))
            termsAll <- unique(termsAll)
            colores <- rep("slategray1", length(termsAll)) # Por defecto, si no hay info de expresion
            names(colores) <- termsAll
            if(!is.null(termExprs))
            {
                colores[which(!is.na(termExprs))] <- "yellow" # Default, if geneExpr available
                colores[which(termExprs>=0.75)] <- "darkgoldenrod1"
                colores[which(termExprs==1)] <- "#FF1010FF" # "darkred"
                colores[which(termExprs<=-0.75)] <- "darkolivegreen1"
                colores[which(termExprs==-1)] <- "#10CF10FF" # green
            }
            
            # Identify terms from each ontology
            gos<- list()
            gos[["BP"]] <- colores[names(colores) %in% names(parentsDB[["BP"]])]
            gos[["MF"]] <- colores[names(colores) %in% names(parentsDB[["MF"]])]
            gos[["CC"]] <- colores[names(colores) %in% names(parentsDB[["CC"]])]
            
            # Plot
            fileNames[[cl]] <- list()
            for(db in c("BP","MF","CC")) 
            {
                if(length(gos[[db]])>0)
                {
                    fileName <- paste(paste("GO_",grPrefix, cl,"_",db,".png", sep=""),sep="")                    
                    tColor <- NULL
                    if(!is.null(termExprs)) tColor <- gos[[db]]
                    leaves[[cl]] <- c(leaves[[cl]], plotGoAncestors(names(gos[[db]]), tColor=tColor, ontology=parentsDB[[db]], fileName=fileName, height=1000, returnLeaves=TRUE)[,1])
                    fileNames[[cl]][[db]] <- paste(folder,fileName, sep="")
                }
            }
        }
        
        names(fileNames) <- paste(clusterColumn, names(fileNames))    
        ret <- list(fileNames=fileNames,leaves=leaves)
        return(ret)
    }else{
        return("") # No GO terms
    }
}

# Main function
createHtml <- function(htmlFileName, feaResults, jobName, tablesGenes, tablesTerms, tool, queryArgs, filterAttribute, filterOperator, filterThreshold, geneExpr=NULL, plotExpression=NULL, onlyGoLeaves=TRUE, plotGoTree=TRUE, plotKeggPw=TRUE)
{
    #####################################################################################################
    ####################################   Initializations   ############################################

    data("groupTypes", envir = environment())
    groupTypes<- get("groupTypes", envir  = environment())
    grType <- names(feaResults)[which(names(feaResults) %in% c("metagroups", "clusters"))]
        
    # Group types...
    sortMG <- TRUE
    if(is.na(filterAttribute)) sortMG <- FALSE
    rawMetagroups <- feaResults[[grType]]
    if(tool=="gage") 
    {
        essentialRawMetagroups <- rawMetagroups
        rawMetagroups <- readGeneTermSets(feaResults$fileName, tool="gage", simplifyGage=FALSE)$clusters
    }
    
    grType <- rownames(groupTypes)[which(groupTypes[,"feaResultsName"]==grType)]
    
    # GTLinker or other?
    if(is.null(rawMetagroups))    
    {
        rawMetagroups <- feaResults$geneTermSets
        rawMetagroups <- data.frame(cbind(Gtset=1:nrow(rawMetagroups),rawMetagroups))
        colsAvg <- NULL
        if(tool=="topGO") colsAvg <- c("Annotated", "classic")
        rawMetagroups <- clustersTable(rawMetagroups, clusterColumn="Gtset", addKeyWordsTerm=FALSE, colsAvg=colsAvg)
        sortMG <- FALSE
            
        grType <- "Gene-term set"
    }
    if(grType == "Metagroup") 
    {
        onlyGoLeaves <- FALSE
    }

    mgIncMatrix <- tablesGenes[[which(names(tablesGenes)==groupTypes[grType,"tablesName"])]]
    termsNumMg <- apply(tablesTerms[[which(names(tablesTerms)==groupTypes[grType,"tablesName"])]],1,sum)
    grPrefix <- groupTypes[grType,"prefix"]
    
    filteredOut <- tablesGenes$filteredOut
    
    genesMetagroups <- unique(unlist(sapply(rawMetagroups[,"Genes"], function(x) strsplit(as.character(x), split=","))))
    nGenesTotal <- length(genesMetagroups)
    nGenesInNotFilteredMg <- nrow(mgIncMatrix)
    nRawMg <- nrow(rawMetagroups)
    
    # Get metagroup colors
    allMgNames <- c(colnames(mgIncMatrix), filteredOut)
    if(all(!is.na(suppressWarnings(as.numeric(allMgNames))))) allMgNames <- as.numeric(allMgNames)
    colores <- setColors(as.character(sort(allMgNames)))[colnames(mgIncMatrix)]
    
    if(sortMG) 
    {
        globalMetagroups <- rawMetagroups[colnames(mgIncMatrix),]
    }else{
        globalMetagroups <- rawMetagroups[rownames(rawMetagroups) %in% colnames(mgIncMatrix),]
    }

    #### Set folder to save images etc: Same as downloaded feaResults    
    folder <- paste(jobName, .Platform$file.sep, sep="")      
    
    # Get metagroup terms
    keggPlots <- NULL   
    geneIDtype <- feaResults$queryArgs$geneIdType
    if(plotKeggPw && !is.null(geneIDtype)) 
    {
        keggExpression <- geneExpr
        if(is.null(keggExpression)) keggExpression <- setNames(rep(NA, length(genesMetagroups)), genesMetagroups)
        
        keggIDs <- getTerms(feaResults, returnValue="KEGG")
        names(keggIDs) <- NULL
        keggIDs <- unlist(keggIDs)
        keggPlots <- tryCatch( 
                    {
                        keggPlots <- plotKegg(keggIDs, keggExpression, geneIDtype=geneIDtype)
                        keggPlots
                    }, error = function(e) 
                    {
                        plotKeggPw <<- FALSE
                    })
        if(all(keggPlots==FALSE)) keggPlots <- NULL
        if(!is.null(keggPlots))   keggPlots <- setNames(paste(folder, keggPlots, sep=""), names(keggPlots))
    }
    mgTerms <- getMGTerms(globalMetagroups, grType, keggPlots=keggPlots)
    
    # Add asterisc to terms in several metagroups:
    termsMultipleMg <- termsNumMg[which(termsNumMg>1)]
    if(length(termsMultipleMg)>0)
    {
        mgTerms <- sapply(mgTerms, function(x) 
        {
            tmp <- rownames(x) %in% names(termsMultipleMg)
            if(any(tmp))  x[tmp,"TermDescription"] <- paste(x[tmp,"TermDescription"], " *", sep="")
            x
        })
    }
    
    # Go trees:
    notFilteredGeneTermSets <- feaResults$geneTermSets[!feaResults$geneTermSets[,1] %in% filteredOut,]
    
    goTrees_global <- NULL
    if(grType=="Gene-term set")  
    {
        # Individual gene-term sets: Plot global GO tree
        goLinks_full <- ""
        goTrees_global <- goTreeLinks(gtSets=cbind("Gene-term set"=1,notFilteredGeneTermSets), grPrefix, geneExpr, folder, plotGoTree)

    }else
    {
        # Clusters or metagroups: plot GO tree by cluster/mg
        goLinks_full <- goTreeLinks(notFilteredGeneTermSets, grPrefix, geneExpr, folder, plotGoTree)
    }
    
    
    leaves <- NULL
    goLinks <- NULL
    if(is.list(goLinks_full))
    {
        if(onlyGoLeaves) leaves <- goLinks_full[["leaves"]]
        goLinks <- goLinks_full$fileNames
    }    
    # Terms table to write in html: 
    boldTerms <- NULL
    if(tool=="gage") boldTerms <- essentialRawMetagroups$Terms
    termsTables <- buildTermsTable(mgTerms, leaves, boldTerms)                                                                        
    
    # Copiar CSS a la carpeta actual...
    cssDir <- file.path(system.file('css', package='FGNet'))
    cssFile <- paste(cssDir, "functionalNetworks.css", sep=.Platform$file.sep)
    file.copy(cssFile, ".")
    
    # Add representative term to DAVID legend
    legendMg <- TRUE
    if(nrow(globalMetagroups)>15) legendMg <- FALSE
    if(legendMg)
    {
        legendMg <- NULL
        if(grType ==  "Cluster") legendMg <- paste(keywordsTerm(mgTerms), ", ...", sep="")
    }
    
    #####################################################################################################
    #################################### Create Plots  ###################################################
    
    iGraphFile <- "iGraph.RData"
    adjMatFile <- "adjMat.RData"
    networkPlot <- "nwFunctionalNetwork.png"
    networkPlot2a <- "nwIntersection_kk.png"
    networkPlot2b <- "nwIntersection_circle.png"
    networkTerms <- "nwTerms.png"
    nwStatsPlot <- "nwStats.png"
    distancePlot <- "plot_Distance.png"
    simplifiedTable <- "simplifiedTermsTable.txt"
    
    # Set node and label size
    numNodes <- sum(dim(mgIncMatrix))
    if(numNodes < 150)
    {
        vSize <- 12
        vLabelCex <- 0.9
    }
    if(numNodes >= 150)
    {
        vSize <- 10
        vLabelCex <- 3/4
    }    
    if(numNodes > 300)
    {
        vSize <- 8
        vLabelCex <- 2/3
    }    
    
    png(networkPlot, width = 800, height = 800)
    retGraph <- functionalNetwork(tablesGenes, plotType="default", plotOutput="static", vSize=vSize, vLabelCex=vLabelCex, legendText=legendMg, geneExpr=geneExpr, plotExpression=plotExpression)
    dev.off()
    png(nwStatsPlot, width = 800, height = 500)
    nwStats <- analyzeNetwork(tablesGenes, fNw=retGraph, plotOutput=TRUE, colors=colores)
    dev.off()
    
    png(networkPlot2a, width = 800, height = 800)
    intersectionGraph <- functionalNetwork(tablesGenes, plotType="bipartite", plotOutput="static", vLayout="kk", vSize=vSize, vLabelCex=vLabelCex, legendPrefix=grPrefix, geneExpr=geneExpr, plotExpression=plotExpression, plotTitle=paste("Genes in several ", tolower(grType), "s", sep=""))
    dev.off()
    if(!is.null(intersectionGraph))
    {
        tablesGenes$iGraph <- c(retGraph$iGraph, intersectionGraph$iGraph)
        retGraph$adjMat <- c(retGraph$adjMat, intersectionGraph$adjMat)
        
        png(networkPlot2b, width = 800, height = 800)
        functionalNetwork(tablesGenes, plotType="bipartite", plotOutput="static", vLayout="circle", vSize=vSize, vLabelCex=vLabelCex, legendPrefix=grPrefix, geneExpr=geneExpr, plotExpression=plotExpression,plotTitle=paste("Genes in several ", tolower(grType), "s", sep=""))
        dev.off()
    }
    
    if(length(termsMultipleMg)>0)
    {
        png(networkTerms, width = 800, height = 800)
        functionalNetwork(tablesTerms, plotType="bipartite", plotAllMg=FALSE, vSize=vSize, vLabelCex=vLabelCex, legendPrefix=grPrefix, plotTitle=paste("Terms in several ", tolower(grType), "s", sep="")) 
        dev.off()
    }
    
    png(distancePlot, width = 600, height = 600)
    distanceMatrix <- clustersDistance(tablesGenes, colores)
    dev.off()

    iGraph<-retGraph$iGraph
    save(iGraph, file=iGraphFile)
    adjMat<-retGraph$adjMat
    save(adjMat, file=adjMatFile)
    
    #####################################################################################################
    #################################### Create HTML  ###################################################
    p=openPage(htmlFileName, link.css=paste(folder, "functionalNetworks.css", sep=""))
    
    # Header
    hwrite('Functional Gene Networks', p, heading=1)
    
     # Query parameters
    createHTML_FEA_header(p, tool, queryArgs, feaResults, folder)
    
    # Results header                
    hwrite("Results: ", p, heading=2)
    
    #hwrite("Basic stats: ", p, br=TRUE)
    hwrite(paste("Number of ",tolower(grType),"s: ",sep=""), p, class='InfoLabel', br=FALSE)
    hwrite(nRawMg, p, br=TRUE)
    hwrite(paste("Number of genes included in all ",tolower(grType),"s: ",sep=""), p, class='InfoLabel', br=FALSE)
    hwrite(nGenesTotal, p, br=TRUE)
    hwrite(paste("Number of genes included in non-filtered ",tolower(grType),"s: ",sep=""), p, class='InfoLabel', br=FALSE)
    hwrite(nGenesInNotFilteredMg, p, br=TRUE)
    hwrite(paste("Filtered ",tolower(grType),"s ",sep=""), p, class='InfoLabel', br=FALSE)
    if(!is.na(filterAttribute)) hwrite(paste(" (", filterAttribute, " ",filterOperator, " ", filterThreshold, ")", sep=""), p, br=FALSE)
    hwrite(": ", p, br=FALSE)
    if(length(filteredOut)>0)  hwrite(paste(grPrefix, filteredOut, sep="", collapse=", "), p, br=TRUE) 
    else hwrite("None", p, br=TRUE) 

    hwrite("<br/>Exports: ", p, br=TRUE)
    hwrite("Functional network: ", p)
    #if(!is.null(pdfName)) hwrite('[PDF] ', p, link=paste(pdfName, sep="" ), br=FALSE)
    hwrite("iGraph (.RData)", p, link=paste(folder,iGraphFile, sep="" ), br=FALSE)
    hwrite(", ", p)
    hwrite("adjacency matrices (.RData)", p, link=paste(folder,adjMatFile, sep="" ), br=TRUE)
    
    if(grType != "Gene-term set")
    {
        hwrite(paste("Simplified ", tolower(grType),"-terms table (shown below): ",sep=""), p)
        hwrite(simplifiedTable, p, link=paste(folder,simplifiedTable, sep="" ), br=FALSE)
    }
    
    # Functional Network plot
    hwrite("Functional Network: ", p, heading=3)
    subImages <- ""
    closeButton <- '<a href="#close" title="Close" class="close">X</a>'
    if(!is.null(intersectionGraph))
    {        
        subImages <- c(hwrite(paste('Genes in several ', tolower(grType),'s: ', sep=""), class='InfoLabel', br=TRUE),
                         hwrite('(Different layouts)'),
                         
                         hwriteImage(paste(folder,networkPlot2a,sep=""), class='intersectionNw', link='#nwIntersection_kk', br=FALSE),
                         hwrite(paste('<div id="nwIntersection_kk" class="modalDialog"><div>',closeButton, '<img border="0" src="',folder,networkPlot2a,'" width="100%"></div></div>', sep="")),
                                                              
                         hwriteImage(paste(folder,networkPlot2b,sep=""), class='intersectionNw', link='#nwIntersection_circle', br=FALSE),
                         hwrite(paste('<div id="nwIntersection_circle" class="modalDialog"><div>',closeButton, '<img border="0" src="',folder,networkPlot2b,'" width="100%"></div></div>', sep="")))

        
        if(length(termsMultipleMg)>0) 
        {
            subImages <- c(subImages,
                         hwrite(paste('Terms in several ', tolower(grType),'s: ', sep=""), class='InfoLabel', br=TRUE),
                         hwrite('(Marked with * in the table)'),
                         hwriteImage(paste(folder,networkTerms,sep=""), class='intersectionNw', link='#nwTerms', br=FALSE),
                         hwrite(paste('<div id="nwTerms" class="modalDialog"><div>',closeButton, '<img border="0" src="',folder,networkTerms,'" width="100%"></div></div>', sep="")))
        }
    
    }
    imageTable <- c(hwriteImage(paste(folder,networkPlot,sep=""), class='network', link='#functionalNetwork', br=FALSE),
                                    hwrite(paste('<div id="functionalNetwork" class="modalDialog"><div>',closeButton, '<img border="0" src="',folder,networkPlot,'" width="100%"></div></div>', sep="")),
                                    hwrite(subImages, border=0, dim=c(length(subImages),1), class="ImageTable"))
    hwrite(c(imageTable), p, border=0, class="ImageTable")
    
    clustSectionHeader <- ""
    if(sortMG) clustSectionHeader <- paste(" (sorted by ",filterAttribute,")", sep="")
    clustSectionHeader <- paste(grType, "s", clustSectionHeader,": ", sep="")
    
    hwrite(clustSectionHeader, p, heading=3, br=FALSE)
    showHideTerms <- '<span id="switchTerms" class="switcher switchTerms"></span><span id="contentTerms">'
    hwrite(showHideTerms, p, br=TRUE)
        
    # Terms table
    termsTable <- NULL
    txtTermsTable <- NULL
    goLinks <- goLinks[names(termsTables$termLinks)]
    for(mg in 1:length(termsTables$termLinks))
    {
        mgName <- paste("<b>",names(termsTables$termLinks)[mg],"</b>")
        attrs <- NULL
        if(tool=="GeneTerm Linker")
        {
            attrs <- c(paste("Silhouette: ", round(globalMetagroups[mg, "Silhouette Width"], 2), sep=""),
                                 paste("P-value: ", signif(globalMetagroups[mg, "pValue"],2),sep=""))        
        }        
        if(tool=="DAVID")
        {    
            attrs <- c("", # "Silhouette  only in GtLinker
                                 paste("Score: ", signif(as.numeric(as.character(globalMetagroups[mg, "EnrichmentScore"])),2),sep=""),
                                 paste("minPval: ", signif(as.numeric(as.character(globalMetagroups[mg, "minPval"])),2),sep=""))
        }    
        if(tool=="topGO")
        {    
            attrs <- c("", # "Silhouette  only in GtLinker
                       paste("Num annotated: ", signif(as.numeric(as.character(globalMetagroups[mg, "Annotated"])),2),sep=""),
                       paste("P-val: ", signif(as.numeric(as.character(globalMetagroups[mg, "classic"])),2),sep=""))
        }    
        if(tool=="gage")
        {    
#            attrs <- c("", # "Silhouette  only in GtLinker
#                        paste("Score: ", signif(as.numeric(as.character(globalMetagroups[mg, "EnrichmentScore"])),2),sep=""),
#                        paste("minPval: ", signif(as.numeric(as.character(globalMetagroups[mg, "minPval"])),2),sep=""))
        }    
        
        
        
        
        tmpGenes <- sort(unlist(strsplit(as.character(globalMetagroups[mg, "Genes"]), ",")))
        if(!is.null(geneExpr))
        {
            ngUP <- sum(geneExpr[tmpGenes]>0, na.rm=TRUE)
            ngDW <- sum(geneExpr[tmpGenes]<0, na.rm=TRUE)
        }
        nGenes <- length(tmpGenes)
        tmpGenes <- paste(tmpGenes, collapse=", ")
        linkGenes <- paste("Genes: ", '<a title="', tmpGenes,'" class="tooltip"><span title="Genes">', nGenes, '</span></a>',sep="")             #href="#close"
        if(!is.null(geneExpr)) 
        {
                linkGenes <- paste(linkGenes, " (", sep="")
                if(ngUP>0) linkGenes <- paste(linkGenes, ngUP, "up", sep="")
                if(ngUP>0 && ngDW>0) linkGenes <- paste(linkGenes, ",", sep="")
                if(ngDW>0) linkGenes <- paste(linkGenes, ngDW, "dw", sep="")
                linkGenes <- paste(linkGenes, ")", sep="")
        }
        txtGenes <- paste(nGenes, " genes: ", tmpGenes, sep="")
        attrsTxt <- c(attrs, txtGenes)
        attrs <- c(attrs, linkGenes)                     
        
        treesLink <- ""
        if(!is.null(goLinks[[mg]]) && all(!is.na(goLinks[[mg]])))
        {
            treesLink <- paste("<div id='thumbnail'>") 
            for(db in names(goLinks[[mg]]))
            {
                if(!is.null(goLinks[[mg]][[db]]))    treesLink <- paste(treesLink, paste("<a href='",goLinks[[mg]][[db]],"'>[",db, "]<img src='",goLinks[[mg]][[db]],"'/></a>", sep=""), collapse=" ")        
            }
            treesLink <- paste(treesLink, "</div>")
            treesLink <- cbind("<div style='text-align:right'>GO terms tree: </div>",treesLink)
        }
        terms <- rbind(termsTables$termLinks[[mg]], treesLink)
        dimnames(terms) <- NULL
        
        termsTable <- c(termsTable, hwrite(c(mgName, attrs), class='mgAttr', table=TRUE, border=0), 
                                        hwrite(terms, table=TRUE, col.class=list("termsRow", "annotRow")))
        txtTermsTable <- rbind(txtTermsTable, paste(c(names(termsTables$noLinks)[mg], attrsTxt), collapse="\t"), rbind(termsTables$noLinks[[mg]]))
    }
    
    #bgCols <- c(strsplit(paste(substr(colores, 1, 7)[rownames(globalMetagroups)], "#FFFFFF", collapse=" "), split=" ")[[1]],"#FFFFFF")
    borderCols <- paste('border-color:', strsplit(paste(substr(colores, 1, 7), substr(colores, 1, 7), collapse=" "), split=" ")[[1]], sep="")
    hwrite(termsTable, p, border=1, class=rep(c('mgHeader', 'Terms'), length(termsTable)/2), row.style=borderCols, dim=c(length(termsTable),1))
    
    # Explanation of asterisc (there are terms in several metagroups):
    if(length(termsMultipleMg)>0) 
    {
        hwrite(paste("<i>* Terms marked with an asterisc are in several ", tolower(grType),"s.</i>", sep=""), p, br=TRUE)
        #hwriteImage(networkTerms, p, link=c(networkTerms), br=TRUE,border=0)    
    } else {
        hwrite(paste("<i>(There are no terms in more than one ", tolower(grType),")</i>", sep=""), p, br=TRUE)
    }

    if(!is.null(boldTerms) && (tool=="gage")) hwrite(paste("Bold indicates <i>essential</i> term.", sep=""), p, br=TRUE)


    hwrite('</span>', p, br=TRUE) # Section end (Terms table)
    
    
    if(!is.null(goTrees_global) && all(goTrees_global!=""))
    {
        
        goTreesImages <- sapply(goTrees_global$fileNames[[1]], function(image) hwriteImage(image, link=image, br=FALSE, class="GoTreesGlobal"))
        
        hwrite(paste('Global GO trees:', sep=""), p, class='InfoLabel', br=TRUE)
        hwrite(goTreesImages, p, border=0, dim=c(1,length(goTreesImages)), class="ImageTable")    
    }

    hwrite("Network analysis: ", p, heading=3, br=FALSE)
    hwrite('<span id="switchNwAnalysis" class="switcher switchNwAnalysis"></span><span id="contentNwAnalysis">', p, br=FALSE)
    
    hwriteImage(paste(folder,nwStatsPlot,sep=""), p, class='nwStats', link=paste(folder,nwStatsPlot,sep=""), br=TRUE)
    
    # Clustering coeficient. 
    hwrite("<b>Clustering coeficient / Transitivity:</b> ", p, br=TRUE)
    hwrite('<div class="transitivity">', p, br=FALSE)
    hwrite(as.matrix(round(nwStats$transitivity, digits=3)), p, table=TRUE, border=0, br=TRUE)
    hwrite('</div>', p, br=FALSE)  
    # hwrite("<i>(Transitivity measures the probability that the adjacent vertices of a vertex are connected)</i>", p, br=TRUE)

    hwrite("<br><b>Inter-modular hubs (whole network):</b> ", p, br=TRUE)
    hwrite(paste(nwStats$hubsList$Global, collapse=", "), p, br=TRUE)
    
    hwrite(paste("<br><b>Intra-modular hubs (within each ", tolower(grType),"):</b>", sep=""), p, br=TRUE)
    intraModularHubs <- nwStats$hubsList[names(nwStats$hubsList)!="Global"]
    # Lista
    # hwrite("Hub list:", p, br=TRUE)
    sapply(names(intraModularHubs), function(cl)
    {
        if(length(intraModularHubs[[cl]])>0) hwrite(paste(grType, " ", cl, ": ",paste(intraModularHubs[[cl]], collapse=", "), sep=""), p, br=TRUE)
    })

    # Comunes
    if(length(nwStats$intraHubsCount)>0) 
    {
        hwrite("<br>Most common hubs: ", p, br=FALSE)
        hwrite(names(nwStats$intraHubsCount), p,border=0, class='matrix', br=TRUE)
    }else{
        hwrite("(There are no hubs in more than one cluster)", p, br=TRUE)
    }
    #hwrite(nwStats$hubs, p, table=TRUE, border=0, class='matrix', br=TRUE)
    hwrite('</span>', p, br=TRUE)  # Section end (contentNwAnalysis)
    
    # Distances Plot    
    if(!is.null(distanceMatrix))
    {   
        hwrite(paste("Distances between ", grType,"s: ", sep=""), p, heading=3, br=FALSE)
        hwrite('<span id="switchTerms" class="switcher switchTerms"></span><span id="contentTerms">', p, br=FALSE)
        
        colnames(distanceMatrix) <- paste(grPrefix, colnames(distanceMatrix), sep="")
        rownames(distanceMatrix) <- paste(grPrefix, rownames(distanceMatrix), sep="")
        distanceMatrix <- round(distanceMatrix,2)
        
        if(ncol(distanceMatrix) < 25) 
        {
            subImages <- c(hwrite('Distance matrix: ', class='InfoLabel', br=TRUE),
                                     hwrite(distanceMatrix, table=TRUE, border=0, class='matrix', br=FALSE))
        }else{
            subImages <- ""
        }
        imageTable <- c(hwriteImage(paste(folder,distancePlot,sep=""), class='distance', link=paste(folder,distancePlot,sep=""), br=TRUE), hwrite(subImages, border=0, dim=c(length(subImages),1)))
        
        hwrite(c(imageTable), p, table=TRUE, border=0) #, class="ImageTable"
        hwrite('</span>', p, br=TRUE)  # Section end (Distances)
    }
    closePage(p)
    
    write.table(txtTermsTable, file=simplifiedTable, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    browseURL(htmlFileName) # Open html
}
