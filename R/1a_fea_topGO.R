
# TOP GO
######################################################
# For help/vignette:
# genesUniverse <- refList(refPackage, geneIdType)
# geneID2GO <- buildDatabases(orgPackage, geneIdType)$geneID2GO  # To execute multiple queries...
# http://www.bioconductor.org/packages/release/BiocViews.html#___ChipDb
# 
# # Arguments
# geneList <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
#               "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
#               "GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
#               "LSM6", "LSM7", "LSM8", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
#               "PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3",
#               "RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73",
#               "SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")
# annotations=c("GO_BP","GO_MF","GO_CC")
# genesUniverse <- NULL
# organism="Sc"
# geneIdType <- "GENENAME"
# nodeSize <- 5  # 1 (no prune), more stable: 5-10
# pValThr <- 0.01
# 
# geneID2GO=NULL
# refPackage=NULL
# nodeSize=5
# pValThr=0.01
# testStat=NULL
# jobName=NULL
# 
# results <- fea_topGO(geneList, geneIdType=geneIdType, organism="Sc", annotations=annotations, refPackage=refPackage, nodeSize=5, pValThr=0.01, testStat=NULL)



fea_topGO <- function(geneList, geneIdType="ENSEMBL", geneLabels=NULL, organism="Hs", annotations=c("GO_BP","GO_MF","GO_CC"), genesUniverse=NULL, refPackage=NULL, geneID2GO=NULL, nodeSize=5, pValThr=0.01, testStat=NULL, jobName=NULL)
{
    if(!loadInstPkg("topGO")) stop("Package topGO is not available.")
    if(!loadInstPkg("GO.db")) stop("Package GO.db is for fea_topGO.")
    groupGOTerms(where=.GlobalEnv)

    ###############################
    # Check arguments 
    if(!is.null(nodeSize) && !is.numeric(nodeSize)) stop("nodeSize should be numeric.")
    if(!is.null(pValThr) && !is.numeric(pValThr)) stop("pValThr should be numeric.")
    if(!is.null(jobName) && !is.character(jobName)) stop("jobName should be character.")
    
    ###############################
    # Prepare jobName / folder
    if(is.null(jobName)) jobName <- paste(sample(100000:999999,size=1), "_topGO",sep="")
    folder <- jobName
    
    # Create folder
    if((!file.exists(file.path(folder))))
    {
        dir.create(file.path(folder))        
    }
    currWD <- getwd()
    setwd(folder)
    
    ###############################
    # Get organism package
    data("organisms", envir = environment())
    organisms<- get("organisms", envir  = environment())
    
    if(!is.character(organism)) stop("Organism should be the name of an organism db package or one of the rownames in the following table: data(organisms)")
    if(organism %in% rownames(organisms)) 
    {
        orgPackage <- organisms[organism,"orgPackage"]
    } else 
    {
        orgPackage <- organism
    }

    # Get background list (genes universe)
    if(!is.null(refPackage) && !is.null(genesUniverse)) stop("Please provide either 'genesUniverse' or 'refPackage'.")
    if(is.null(refPackage)) refPackage <- orgPackage
    if(is.null(genesUniverse)) genesUniverse <- refList(refPackage, geneIdType)
    
    # Get GO gene sets
    if(!is.null(orgPackage) && !is.null(geneID2GO)) stop("Please provide either 'orgPackage' or 'geneID2GO'.")
    if(any(!annotations %in% c("GO_BP","GO_MF","GO_CC"))) stop('"annotations" should be "GO_BP", "GO_MF" and/or "GO_CC"')
    if(is.null(geneID2GO)) geneID2GO <- buildGeneSets(orgPackage, geneIdType, annotations=annotations)$genes2terms    

    
    # Reshape input list
    recogPercent <- sum(geneList %in% genesUniverse)/length(geneList)
    if(recogPercent<0.9) 
    {
        if(recogPercent < 0.1) stop("There are no genes from the gene list in the genes universe.\nMake sure the arguments 'genesUniverse', 'refPackage', 'organism' and 'geneIdType' are correct.")
        warning(paste("Only ",recogPercent*100, "% of gene ids were recognized.",sep=""))
    }
    inputGeneList <- factor(as.integer(genesUniverse %in% geneList))
    names(inputGeneList) <- genesUniverse
    
    ## topGo
    annotations <- gsub("GO_", "", annotations)
    names(geneID2GO) <- gsub("GO_", "", names(geneID2GO))
    
    GOdata <- list()
    goResTable <- list()
    if(is.null(testStat)) testStat <- new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
    for(ont in annotations)
    {
        GOdata[[ont]] <- new("topGOdata", ontology=ont, allGenes=inputGeneList, annot=annFUN.gene2GO, gene2GO=geneID2GO[[ont]], nodeSize=nodeSize)
        
        # Get significant terms and build table: 
        resultFisher <- getSigGroups(GOdata[[ont]], testStat)
        goResTable[[ont]] <- GenTable(GOdata[[ont]], classic=resultFisher, ranksOf="classic", topNodes=sum(resultFisher@score<pValThr), numChar=100)
        
        # Check size?
        # Get genes in each term:
        goResTable[[ont]] <- cbind(Ont=ont, goResTable[[ont]], Genes=sapply(goResTable[[ont]]$GO.ID,function(go) 
        {
            gTerm <- genesInTerm(GOdata[[ont]], go)[[1]]
            paste(gTerm[gTerm %in%  sigGenes(GOdata[[ont]])], collapse=",")
        }))    
    }
    geneTermSets <- do.call(rbind,goResTable)
    
    ###############################
    # Adaptar formato
    fileName <- paste(jobName, ".txt", sep="")
    
    
    # Equivalent to format_results(goResTable, newFileName=fileName, clusterCol=NULL, geneCol=NULL, geneSep=NULL, termDescCol="Term", termIDCol="GO.ID", termCatCol=NULL, termCat="GO", termSep=NULL), but deltes extra columns
    colnames(geneTermSets)[which(colnames(geneTermSets)=="Term")] <- "Terms"
    geneTermSets$Terms <- paste(geneTermSets$GO.ID, capitalize(geneTermSets$Terms), sep=":")
    geneTermSets <- geneTermSets[,!colnames(geneTermSets) %in% "GO.ID"] 
    
    # Replace Gene ID?
    geneTermSets <- addGeneLabels(geneTermSets, geneLabels)

    # Save file: 
    ret <- format_results(geneTermSets, newFileName=fileName, tool="topGO")
    queryArgs <- list(queryArgsAsCharacter(match.call()))
    setwd(currWD)
    invisible(c(queryArgs=queryArgs,ret))
}
