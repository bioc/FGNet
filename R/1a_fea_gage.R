# ##### Arguments
# library(gage); data(gse16873)
# eset <- gse16873
# refSamples <- grep('HN',colnames(gse16873), ignore.case =T)
# compSamples <- grep('DCIS',colnames(gse16873), ignore.case =T)
# 
# geneIdType <- "ENTREZID"
# annotations <- "REACTOME" # c("GO_BP","GO_MF","GO_CC","KEGG","REACTOME")   , Set to NULL to use the geneSets as is (not named with annotation)
# organism <- "Hs" 
# # required if no geneSets are provided 
# geneSets <- NULL #entrezGS$terms2genes$REACTOME    # Named by annotation (optional, to filter)
#   library(gage)
#   C2geneSets <- readList("c2.cp.v4.0.symbols.gmt")
# geneSets <- lapply(geneSets[c(1,5)], function(x) sample(x,100)) 

# jobName <-  NULL
# 
# sameDirection <- FALSE    
# compareType <- "as.group"
# # 'as.group', group-on-group comparison between ref and samp; 
# # 'unpaired' (used to be '1on1'), one-on-one comparison between all possible ref and samp combinations, although the original experimental design may not be one-on-one paired; 
# # '1ongroup', comparison between one samp column at a time vs the average of all ref columns.
# # ... More arguments to pass to gage




##############################################
fea_gage <- function(eset, refSamples, compSamples, geneIdType, geneLabels=NULL, organism="Hs", annotations=c("GO_BP","GO_MF","GO_CC","REACTOME"), geneSets=NULL, sameDirection=FALSE, onlyEssentialTerms=TRUE, compareType="as.group", jobName=NULL, ...)
{        
    if(!loadInstPkg("gage")) stop("Package gage is not available.")
    
    queryArgs <- list(queryArgsAsCharacter(match.call()))
    
    ######################
    # Check arguments 
    exprsEset <- eset
    if(!is.matrix(eset)) exprsEset <- exprs(eset)
    if(!is.null(jobName) && !is.character(jobName)) stop("jobName should be character.")
    if(any(refSamples %in% compSamples)) stop("refSamples and compSamples should not overlap.")
    if(is.character(refSamples)) refSamples <- which(colnames(eset) %in% refSamples)
    if(is.character(compSamples)) compSamples <- which(colnames(eset) %in% compSamples)
    
    data("organisms", envir = environment())
    organisms<- get("organisms", envir  = environment())
    if(!is.character(organism)) stop("Organism should be the name of an organism db package or one of the rownames in the following table: data(organisms)")
    if(organism %in% rownames(organisms)) 
    {
        orgPackage <- organisms[organism,"orgPackage"]
    } else 
    {
        orgPackage <- organisms
    }
        
    ###############################
    # Prepare jobName / folder
    if(is.null(jobName)) jobName <- paste(sample(100000:999999,size=1), "_gage",sep="")
    folder <- jobName
    
    # Create folder
    if((!file.exists(file.path(folder))))
    {
        dir.create(file.path(folder))        
    }
    currWD <- getwd()
    setwd(folder)
    
    ######################
    # Gene Sets
    # If not provided... calculate
    if(is.null(geneSets))
    {
        if(is.null(organism) || is.null(geneIdType) || is.null(annotations)) stop("Please, provide either 'geneSets' OR the annotations, organism and geneIdType to build new gene-sets.")
        geneSets <- buildGeneSets(organismDb=orgPackage,geneIDtype=geneIdType, annotations=annotations, termLabel="TERM")$terms2genes
        fileName <- "geneSets.RData"
        save(geneSets, file=fileName) 
        message(paste("Gene-sets saved as ", fileName, sep=""))
    }
    
    # Filter the geneSets by annotation
    if(!is.null(annotations) && is.list(geneSets[[1]]))
    {
        if(any(!annotations %in% names(geneSets))) stop("The requested annotations do not match names(geneSets).")
        geneSets <- geneSets[annotations]
        
        # Unlist...
        # Paste annotation into id:
        names(geneSets)[which(names(geneSets)=="REACTOME")] <- "REACT"
        for(annot in names(geneSets)[!grepl("GO",names(geneSets))]) # Non-GO annotations
        {
            names(geneSets[[annot]]) <- paste(annot, ":",names(geneSets[[annot]]),sep="")
        }
        # Unlist
        names(geneSets) <- NULL
        geneSets <- unlist(geneSets, recursive=FALSE)
    }    

    #####################
    # GAGE
    resultGage <- gage(exprsEset, ref=refSamples, samp=compSamples, compare=compareType, gsets=geneSets, same.dir=sameDirection)#, ...)

    fileName <- paste(jobName, "_rawGageResults.RData", sep="")
    save(resultGage, file=fileName) 
    message(paste("Raw GAGE results saved as ", fileName, sep=""))
    
    # GS redundancy
    gageGrouped <- list()
    if(sameDirection)
    {
        gageGrouped[["up"]] <- esset.grp(resultGage$greater, exprsEset, gsets=geneSets, ref=refSamples, samp=compSamples, compare=compareType, same.dir=sameDirection, test4up=TRUE)#, ...)
        gageGrouped[["dw"]] <- esset.grp(resultGage$less, exprsEset, gsets=geneSets, ref=refSamples, samp=compSamples, compare=compareType, same.dir=sameDirection, test4up=FALSE)#, ...)
    }else{
        gageGrouped[["both"]] <- esset.grp(resultGage$greater, exprsEset, gsets=geneSets, ref=refSamples, samp=compSamples, compare=compareType, same.dir=sameDirection)#, ...) # test4up not needed
    }
              
    #####################################
    # Generate geneTermSets table (Convert to FGNet format)
    # sapply(gageGrouped$both, length)
    # essentialSets  setGroups    allSets connectedComponent      overlapCounts       overlapPvals       coreGeneSets 
    #            19        19          65                 19               4225               4225                 65 
    
    # Select the significant terms (gageGrouped[[1]]$allSets) from resultGage (keep the 5 columns with info)
    # [[1]] might be "up" or "both"
    geneTermSets <- NULL  # in case only down
    if(length(gageGrouped[[1]]$allSets) > 0) 
    {
        geneTermSets <- resultGage$greater[rownames(resultGage$greater)%in%gageGrouped[[1]]$allSets,1:5, drop=FALSE]
        geneTermSets <- cbind(dir=names(gageGrouped)[1], geneTermSets) # Add dir
    }
    # Add down:
    if(length(gageGrouped)>1)  # = if(sameDir) 
    {
        geneTermSets <- rbind(geneTermSets, cbind(dir="Down",resultGage$less[rownames(resultGage$less)%in%gageGrouped[["dw"]]$allSets,1:5, drop=FALSE]))
    }    
    
    # Add "Terms"  and "essentialSet" columns
    geneTermSets <- data.frame(essentialSet=rownames(geneTermSets) %in%  unlist(sapply(gageGrouped, function(x) x$essentialSets)),
                               geneTermSets, Terms=rownames(geneTermSets))
    # Add Cluster and Genes...
    tmpGenes <- setNames(rep("", nrow(geneTermSets)), rownames(geneTermSets))
    tmpCluster <- setNames(rep("", nrow(geneTermSets)), rownames(geneTermSets))
    for(updw in names(gageGrouped))
    {
        # Term - Cluster
        tmp <- unlist(mapply(rep, 1:length(gageGrouped[[updw]]$setGroups), sapply(gageGrouped[[updw]]$setGroups, length))) # equivalent to melt(gageGrouped[[updw]]$setGroups)
        tmp <- cbind(Terms=unlist(gageGrouped[[updw]]$setGroups), Cluster=tmp)
        # Add direction to cluster name
        tmp[,"Cluster"] <- paste(updw, as.numeric(tmp[,"Cluster"]), sep="_")
    
        # Add genes
        tmp <- cbind(tmp, Genes=sapply(as.character(tmp[,"Terms"]), function(x)
        {
            genesGtset <- gageGrouped[[updw]]$coreGeneSets[[x]]  ### Genes en lista original??
            paste(as.character(genesGtset), collapse=",")
        }))
        rownames(tmp) <- NULL
        
        tmpCluster[as.character(tmp[,"Terms"])] <- tmp[,"Cluster"]
        tmpGenes[as.character(tmp[,"Terms"])] <- tmp[,"Genes"]
    }
    geneTermSets <- data.frame(Cluster=tmpCluster, geneTermSets, Genes=tmpGenes)
    rownames(geneTermSets) <- NULL
    
    # Sort by cluster
    geneTermSets <- geneTermSets[order(as.numeric(do.call(rbind,strsplit(as.character(geneTermSets$Cluster),"_"))[,2])),,drop=FALSE]
    
    # Replace Gene ID?
    geneTermSets <- addGeneLabels(geneTermSets, geneLabels)
    
    # Write file
    fileName <- paste(jobName, ".txt", sep="")
    write.table(geneTermSets, file=fileName, quote=FALSE, row.names=FALSE, sep="\t")
    message(paste("Gene-term sets table saved as ", fileName, sep=""))
    
    #####################################
    # Calculate gene's FC (for plotting)
    genesFC <- log2(apply(exprsEset[,compSamples],1, mean)/apply(exprsEset[,refSamples],1, mean))
    if(!is.null(geneLabels))
    {
        subNames <- which(names(genesFC) %in% names(geneLabels))
        names(genesFC)[subNames] <- geneLabels[names(genesFC)[subNames]]
    }

    #####################################
    # Load results formatted (cluster table... etc)
    ret <- readGeneTermSets(fileName, tool="gage", simplifyGage=onlyEssentialTerms)
    setwd(currWD)
    invisible(c(queryArgs=queryArgs, ret, genesFC=list(genesFC)))
}






# 
# ###################### EXAMPLE
# matr <- toMatrix(geneTermSets)
# matr <- toMatrix(tablaClusters)
# functionalNetwork(matr)
# #matr <- toMatrix(geneTermSets[grepl("up",geneTermSets$Cluster),])
# 
# #####################
# # Expression
# exprsGenes <- apply(exprsEset[,compSamples],1, mean) - apply(exprsEset[,refSamples],1, mean)
# exprsGenes <- log2(apply(exprsEset[,compSamples],1, mean)/apply(exprsEset[,refSamples],1, mean))
# exprsGenes_Subset <- exprsGenes[rownames(matr$clustersMatrix)]
# 
# # Translate geneID 
# gsym <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(matr[[2]]), columns="SYMBOL",keytype=geneIdType)
# gsymTable <-table(gsym[,"SYMBOL"])
# for(y in names(gsymTable[gsymTable>1]))
# {
#     tmp <- which(gsym[,"SYMBOL"]==y)
#     gsym[tmp,"SYMBOL"] <- paste(y, 1:length(tmp), sep=".")
# }
# # matr
# matr[1:2] <- lapply(matr[1:2], function(x) 
# {  
#     rownames(x) <-  sapply(rownames(x), function(y) gsym[gsym[,geneIdType]==y, "SYMBOL"][1])
#     x
# })
# # expr
# names(exprsGenes_Subset) <- sapply(names(exprsGenes_Subset), function(y) gsym[gsym[,geneIdType]==y, "SYMBOL"][1])
# 
# functionalNetwork(matr,geneExpr=exprsGenes_Subset,legendText=F)
