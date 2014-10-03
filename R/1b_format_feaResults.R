
# ... further argumets to pass to "read.csv"
# termSep: Splits multiple terms in "cell"/line
# termCatCol termIDCol: are pasted to term description, if only one term per line
format_results <- function(fileName, newFileName=NULL, clusterCol=NULL, geneCol=NULL, geneSep=NULL, termDescCol=NULL, termIDCol=NULL, termCatCol=NULL, termCat=NULL, termSep=NULL, tool="Imported text file", simplifyGage=TRUE, ...)
{
    rawResults <- fileName
    if(is.character(fileName) && file.exists(fileName)) rawResults <- read.csv(fileName, ...)
    
    # Format genes
    if(!is.null(geneSep)) rawResults[,geneCol] <- sapply(rawResults[,geneCol], function(x) paste(unlist(strsplit(as.character(x),split=geneSep)), collapse=","))
    if(!is.null(geneCol)) colnames(rawResults)[colnames(rawResults)==geneCol] <- "Genes" 
    
    # Format terms
    # Option A) One term per row, add ID or Category
    termInfo <- !all(is.null(termIDCol), is.null(termCatCol), is.null(termCat))
    if(!is.null(termSep) && termInfo) stop("Term ID and term category can only be used if there is only one term per row/gene-term set (termInfo=NULL).")
    if(!is.null(termCatCol) && !is.null(termCat)) stop("Please provide either termCatCol (column name) or termCat (a common value for all terms).")
    if(!is.null(termCat))
    {
        rawResults <- cbind(rawResults, TermCat=rep(termCat, nrow(rawResults)))
        termCatCol <- "TermCat"
    }
    tempTerm <- NULL
    if(!is.null(termCatCol))
    {
        tempTerm <- paste(rawResults[,termCatCol], ":", sep="")  # (avoid GO:GO:0000... )
        tempTerm[grep("GO:", rawResults[,termIDCol], fixed=TRUE)]<-""
    }
    if(!is.null(termIDCol)) tempTerm <- paste(tempTerm, rawResults[,termIDCol], ":", sep="")  
    if(!is.null(termDescCol)) 
    {
        tempTerm <- paste(tempTerm, rawResults[,termDescCol], sep="")
        if("Terms" %in% colnames(rawResults)) colnames(rawResults)[which(colnames(rawResults)=="Terms")] <- "Terms_raw"
        rawResults <- cbind(rawResults, Terms=tempTerm)
    }    
    # Option B) Several terms per row
    if(!is.null(termSep)) 
    {
        rawResults[,termDescCol] <- sapply(rawResults[,termDescCol], function(x) paste(unlist(strsplit(as.character(x),split=termSep)), collapse=","))
        if(!is.null(termDescCol)) colnames(rawResults)[colnames(rawResults)==termDescCol] <- "Terms" 
    }
    
    # Rename clusters
    if(!is.null(clusterCol)) colnames(rawResults)[colnames(rawResults)==clusterCol] <- "Cluster" 

    if(is.null(newFileName)) 
    {
        splitedFN <- strsplit(fileName,".", fixed=TRUE)[[1]]
        fName <- paste(splitedFN[1:(length(splitedFN)-1)], collapse=".")
        fExt <- NULL
        if(length(splitedFN)>1) fExt <- paste(".", splitedFN[length(splitedFN)], sep="")
        newFileName <- paste(fName, "_formatted", fExt, sep="")
    }

    write.table(rawResults, file=newFileName, quote=FALSE, sep="\t", row.names=FALSE)
    
    invisible(readGeneTermSets(newFileName, tool=tool, simplifyGage=simplifyGage))
}
