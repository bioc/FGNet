# fileName <- "goResTable.txt" 
readGeneTermSets <- function(fileName, tool=NULL, simplifyGage=TRUE)
{
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
        
    # Read file
    sepChar <- "\t"
    geneTermSets <- read.table(fileName, sep=sepChar, header=TRUE, quote="")    
    
    # Create clusters table
    colsAvg <- NULL
    sortBy <- NULL
    if(!is.null(tool))
    {
        colsAvg <- sortBy <- FEA_tools[tool,"DefaultFilter"]
        if(is.na(colsAvg)) colsAvg  <- sortBy <- NULL
        
    }
    if(!is.null(tool) && tool=="gage") 
    {
        # gage:
        keepTermsCol <- NULL
        if(simplifyGage) keepTermsCol <- "essentialSet"
        tablaClusters <- clustersTable(geneTermSets, colsAvg="dir", sortBy=NULL, keepTermsCol=keepTermsCol, addKeyWordsTerm=FALSE)
    }else
    {
        # Default, any other case
        tablaClusters <- clustersTable(geneTermSets, colsAvg=colsAvg, sortBy=sortBy, decreasing=TRUE) # keepTermsCol=NULL, addKeyWordsTerm=TRUE, 
    }
    
    
    # Return
    invisible(list(clusters=tablaClusters, geneTermSets=geneTermSets, fileName=fileName))
}
