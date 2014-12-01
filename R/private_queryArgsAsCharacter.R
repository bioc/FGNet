# functCall <- feaResults$call[[1]]

queryArgsAsCharacter <- function(functCall)
{       
    # As list:
    queryArgs <- formals(eval(as.name(functCall[[1]]))) # Default arguments
    queryArgs <- sapply(names(queryArgs),function(x) # Given arguments (raw value, can be a variable name or not)
    {
        tmp <- functCall[[x]]
        if(!is.null(tmp))
        {
            return(tmp)
        }else{
            return(queryArgs[[x]]) # NULL
        }
    })
    queryArgs <- as.list(c(fun=as.character(as.name(functCall[[1]])), queryArgs))
    
    # Eval/format as string:    
    if(("organism" %in% names(queryArgs)) && is.name(queryArgs$organism))
    {
        # Eval and capitalize
        queryArgs$organism <- eval(queryArgs$organism, envir=parent.frame())
        if(!is.null(queryArgs$organism)) queryArgs$organism <- capitalize(tolower(queryArgs$organism))
    }
    if(("geneList" %in% names(queryArgs)) && is.name(queryArgs$geneList))
    {
        # Eval and collapse
        queryArgs$geneList <- paste(eval(queryArgs$geneList, envir=parent.frame()), sep="", collapse=", ")
    }
    if(("argsWS" %in% names(queryArgs)))
    {
        # Eval
        if(is.name(queryArgs$argsWS)) queryArgs$argsWS <- eval(queryArgs$argsWS, envir=parent.frame())
            
        # Collapse (format)
        if(!is.null(names(queryArgs$argsWS)))
        {
            queryArgs$argsWS <- paste(names(queryArgs$argsWS), queryArgs$argsWS, sep="=", collapse=", ")
        }else{
            queryArgs$argsWS <- paste(queryArgs$argsWS, collapse=", ")
        }
    }

    # Other arguments...    
    argumentsLeft <- names(queryArgs)[!names(queryArgs) %in% c("fun","geneList","organism", "argsWS", "...")] 
    argsEval <- c("serverWeb", "geneIdType", "annotations", "email", "jobName", "alreadyDownloaded", "jobID", "refSamples", "compSamples", "sameDirection", "compareType", "nodeSize", "pValThr", "testStat")
            # NOT eval: argsWS, eset, geneSets, genesUniverse, refPackage, geneID2GO
    
    for(arg in argumentsLeft)
    {
        if((is.name(queryArgs[[arg]])) && (arg %in% argsEval)) queryArgs[[arg]] <- eval(queryArgs[[arg]], envir=parent.frame())     
        if(is.null(queryArgs[[arg]])) queryArgs[[arg]] <- "NULL"
        
        queryArgs[[arg]] <- paste(queryArgs[[arg]], sep="", collapse=", ")
    }
    
    return(queryArgs)
}