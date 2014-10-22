# keggIDs <-  unlist(getTerms(results, returnValue="KEGG"))
# # geneExpr <-  sample(c(rnorm(n=10),NA), length(genesYeast), replace=T); names(geneExpr) <- genesYeast
# geneExpr <- keggExpression

plotKegg <- function(keggIDs, geneExpr, geneIDtype="ENSEMBL", colType=c("continuous", "discrete"))
{
    if(!loadInstPkg("KEGGprofile")) stop("Package KEGGprofile is required to plot KEGG pathways.")
    # Functions used: parse_XMLfile, plot_profile
    # Warning: KEGG.db is deprecated. (KEGG.db is used for downloading all pathways and find significant patways. --> Not used here)
    colType <- colType[1]
        
    if(is.null(keggIDs)) return(NULL)
    
    ### Find organism and load data
    keggPrefix <- substring(keggIDs[1],1,3)
    if(keggPrefix=="map") keggPrefix <- "hsa" # TO DO?

    keggIDs <- as.matrix(unique(data.frame(keggIDs)))[,1] # Para que no se pierdan los nombres
    #if(!is.null(names(keggIDs))) names(keggIDs) <- sapply(names(keggIDs), function(x){ tmp <- strsplit(x,":")[[1]]; return(tmp[length(tmp)])})
    
    if(!"organisms" %in% ls()) 
    {
        data("organisms", envir = environment())
        organisms<- get("organisms", envir  = environment())
    }
    orgDb <- organisms[which(organisms[,"keggPrefix"]==keggPrefix),"orgPackage"]    
    if(!loadInstPkg(orgDb)) stop(paste("Package '", orgDb, "' is not available.", sep=""))
    orgDb <- eval(parse(text=orgDb))
    
    # Check provided ID (for compatibility with report and DAVID)
    if(!geneIDtype %in% columns(orgDb)) geneIDtype <- columns(orgDb)[sapply(columns(orgDb), function(x) grepl(x, geneIDtype))]
    if(length(geneIDtype)!=1) stop("geneIDtype not valid")

    numRecogKeys <- sum(names(geneExpr) %in% keys(orgDb, keytype=geneIDtype))
    if(numRecogKeys == 0)
    {
        # Try with gene symbol or labels, in case geneLabels was used
        geneIDtype <- c("SYMBOL", "GENENAME")[which(c("SYMBOL", "GENENAME") %in% columns(orgDb))][1]
        numRecogKeys <- sum(names(geneExpr) %in% keys(orgDb, keytype=geneIDtype))
    }        
        
    fileNames <- NULL
    if(numRecogKeys > 0)
    {    
        ### Select required ID
        # For most species: Entrez
        geneIDrequired <- "ENTREZID" 
        if(orgDb$packageName=="org.At.tair.db") geneIDrequired <- "TAIR"    # Kegg prefix: "ath"
        if(orgDb$packageName %in% c("org.Pf.plasmo.db", "org.Sc.sgd.db")) geneIDrequired <- "ORF"    # Kegg prefix: "pfa","sce" 
        if(!geneIDrequired %in% columns(orgDb)) stop("Gene ID is not available in organism database.")
        
        newIDs <- select(orgDb, keys=names(geneExpr), columns=geneIDrequired, keytype=geneIDtype)
        newIDs <- newIDs[!is.na(newIDs[,geneIDrequired]),]
        
        geneExprNewIDs <- sapply(unique(newIDs[,geneIDrequired]), function(x) mean(geneExpr[unique(newIDs[newIDs[,geneIDrequired]==x, geneIDtype])])) # Mean in case there are several values for the same ID
        names(geneExprNewIDs) <- unique(newIDs[,geneIDrequired])
            
        colsExpr <- rep("skyblue2", length(geneExprNewIDs)) # Default (No expression info)    
        if(!all(is.na(geneExprNewIDs))) 
        {
            if(colType=="continuous" && length(table(geneExprNewIDs))>2) 
            {
                # "#FF1010" rgb(1,0.063,0.063) # darkred
                # "#10CF10" rgb(0.063,0.81,0.063) # green
                colsExpr <- color.scale(geneExprNewIDs,cs1=c(0.063,1,1), cs2=c(0.81,1,0.063), cs3=c(0,0.7,0.063), na.color="skyblue2", xrange=c(min(c(-1,geneExprNewIDs), na.rm=TRUE),max(c(1,geneExprNewIDs), na.rm=TRUE)))
            }else  # if(colType=="discrete") or levels in expression <=2 (with 2 values color.scale only plots black and white)
            {     
                quantsExpr <- max(abs(geneExprNewIDs), na.rm=TRUE)*c(0.1, 0.8)
                colsExpr[which(!is.na(geneExprNewIDs) & geneExprNewIDs>0)] <- "#FF1010FF" # "darkred"
                colsExpr[which(geneExprNewIDs<=quantsExpr[2])] <- "darkgoldenrod1"
                colsExpr[which(geneExprNewIDs<=quantsExpr[1])] <- "yellow" 
                colsExpr[which(geneExprNewIDs<=-quantsExpr[1])] <- "darkolivegreen1"
                colsExpr[which(geneExprNewIDs<=-quantsExpr[2])] <- "#10CF10FF" # green
            }
        }
        names(colsExpr) <- names(geneExprNewIDs)
        #cols <- unique(names(table(sort(colsExpr)))); pie(rep(1,length(cols)),col=cols)
       
        allIds<- keys(orgDb, keytype=geneIDrequired)
        allG <- rep(0, length(allIds))
        names(allG) <- allIds
        geneExprNewIDs <- geneExprNewIDs[!is.na(geneExprNewIDs)]
        allG[names(geneExprNewIDs)] <- geneExprNewIDs
        
        colsAllG <- rep("white", length(allG))
        names(colsAllG) <- names(allG)
        colsAllG[names(colsExpr)] <- colsExpr
        #table(colsAllG)
        
        keggIDs <- setNames(sapply(keggIDs, function(x) gsub(keggPrefix,"",x)), names(keggIDs))
        
        # Pathway names?
        if(!is.null(names(keggIDs))) 
        {
            pathwayNames <- setNames(names(keggIDs), keggIDs)
        }else
        {
            pathwayNames <- setNames(rep("profile", length(keggIDs)), keggIDs) 
            if(loadInstPkg("KEGG.db"))
            {
                pathwayNames <- unlist(AnnotationDbi::as.list(KEGGPATHID2NAME)[keggIDs])
            }        
        }
        
        # Plot
       
        for(pwId in keggIDs)
        { 
            if(!file.exists(paste(keggPrefix,pwId, ".xml",sep=""))) 
            {
                download_KEGGfile(pathway_id=pwId, species = keggPrefix)
                
                # Uses (webLinks?):
                #             http://www.genome.jp/kegg-bin/download?entry=05130&format=kgml
                #             http://www.genome.jp/kegg/pathway/hsa/hsa05130.png
                # REST API usage is recommended. Is the information the same?    (http://www.kegg.jp/kegg/rest/)
                #             http://rest.kegg.jp/get/hsa05130/image
                #             http://rest.kegg.jp/get/hsa05130/kgml       
                
                # download.file(paste("http://rest.kegg.jp/get/", keggPrefix,pwId, "/kgml",sep=""),  paste(keggPrefix,pwId, ".xml",sep=""))
                # download.file(paste("http://rest.kegg.jp/get/", keggPrefix,pwId, "/image",sep=""),  paste(keggPrefix,pwId, ".png",sep=""))            
            }
            XML2database <- parse_XMLfile(pathway_id=pwId, species = keggPrefix)
            try(temp <- plot_profile(cbind(allG), bg_col=cbind(colsAllG),text_col="black", border_col="black", type="bg", pathway_name=paste(keggPrefix,pwId, sep=""), KEGG_database=XML2database, magnify=1, pathway_min=1))    
            file.remove(paste(keggPrefix,pwId, ".xml",sep=""))
            file.remove(paste(keggPrefix,pwId, ".png",sep=""))
            
            pwName <- pathwayNames[pwId]
            fileName <- gsub(" ", "_",gsub("_profile_bg",paste("_",pwName,sep=""),paste(keggPrefix,pwId, "_profile_bg.png",sep="")))
            names(fileName) <- pwId
            file.rename(from=paste(keggPrefix,pwId, "_profile_bg.png",sep=""),to=fileName)
            fileNames <- c(fileNames, fileName)
        }
    }
    invisible(fileNames)
}