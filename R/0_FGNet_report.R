
#ADD? jobName=NULL,  geneLabels=NULL,
FGNet_report <- function (feaResults, geneExpr=NULL, plotExpression="border", onlyGoLeaves=TRUE, plotGoTree=TRUE, filterAttribute=NULL, filterOperator=NULL, filterThreshold=NULL)
{
    #####################################################################################################
    ####################################   Check arguments   ############################################

    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    # Identify tools and extract arguments
    tool <- "Imported text file"
    queryArgs <- NULL
    jobName <- NULL
    if("queryArgs" %in% names(feaResults))
    {
        tool <- rownames(FEA_tools)[which(FEA_tools[,"Function"]==feaResults$queryArgs$fun)]
        queryArgs <- feaResults$queryArgs
        jobName <- queryArgs$jobName
    }
    
    # if(plotKeggPw)
    # {
    #     if(tool=="GeneTerm Linker")
    #     {
    #         plotKeggPw <- FALSE
    #         warning("Local plots of KEGG pathways are not available for Gene-Term Linker. Links to the website will be used instead.")
    #     }
    # }
    
    # Get info:
    if(is.null(filterAttribute)) 
    {
        filterAttribute <- FEA_tools[tool,"DefaultFilter"]
        if(is.null(filterOperator))   filterOperator <- FEA_tools[tool,"DefaultFiltOperator"]
        if(is.null(filterThreshold))  filterThreshold <- FEA_tools[tool,"DefaultFiltThreshold"]
    }
    
    if(is.null(geneExpr) && ("genesFC" %in% names(feaResults))) geneExpr <- feaResults$genesFC
    
    #####################################################################################################
    ####################################   Common to all tools   ######################################## 
    tablesGenes <- fea2incidMat(feaResults, filterAttribute=filterAttribute, filterOperator=filterOperator, filterThreshold=filterThreshold, key="Genes") 
    tablesTerms <- suppressMessages(fea2incidMat(feaResults, filterAttribute=filterAttribute, filterOperator=filterOperator, filterThreshold=filterThreshold, key="Terms"))
    
    #####################################################################################################
    ####################################   Generate  HTML    ############################################
    ###############################
    # Prepare jobName / folder
    # if(is.null(jobName)) jobName <- paste(sample(100000:999999,size=1), "_", tool,sep="")
    if(is.null(jobName) || jobName=="NULL") jobName <- getJobName(feaResults$fileName)    
    
    # Create folder
    folder <- jobName
    if((!file.exists(file.path(folder))))
    {
        dir.create(file.path(folder))        
    }
    currWD <- getwd()
    setwd(folder)
    
    tryCatch({
        htmlFileName <- paste(currWD, .Platform$file.sep, jobName, ".html", sep="") 
        createHtml(htmlFileName=htmlFileName, feaResults=feaResults, jobName=jobName, tablesGenes=tablesGenes, tablesTerms=tablesTerms,  # Data
                   tool=tool, queryArgs=queryArgs, # Query info
                   filterAttribute=filterAttribute, filterOperator=filterOperator, filterThreshold=filterThreshold, # Filter info
                   geneExpr=geneExpr, plotExpression=plotExpression, onlyGoLeaves=onlyGoLeaves, plotGoTree=plotGoTree)  # HTML options
        setwd(currWD)
    }, error = function(e) {
        setwd(currWD)
        stop(e)
    })
}
