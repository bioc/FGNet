createHTML_FEA_header<- function(p, tool, queryArgs, feaResults, folder)
{
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
        
    if(tool!="Imported text file") 
    {
        hwrite(paste("Functional enrichment with ", tool, "",sep=""), p, heading=2)
    }else{
        hwrite("Functional enrichment", p, heading=2)
        hwrite("Report generated from text file with FEA results.", p, heading=2)
    }
    if(!is.null(queryArgs)) hwrite('<span id="switchQueryInfo" class="switcher switchQueryInfo"></span><span id="contentQueryInfo">', p, br=FALSE)
    
    ############################################
    ## Tool
    serverWeb <- FEA_tools[tool,"URL"]
    if(!is.na(serverWeb)) 
    {      
        hwrite("Server: ", p, class='InfoLabel', br=FALSE)
        if("serverWeb" %in% names(queryArgs)) serverWeb <- queryArgs$serverWeb
        hwrite(serverWeb, link=paste(serverWeb, '" target="_blank', sep=""),p, br=TRUE)
    }
       
    ############################################
    ## txt file
    if(tool == "GeneTerm Linker") #### GTLinker !!!!
    {
        if(!is.null(queryArgs$jobID))
        {
            hwrite("Job ID: ", p, class='InfoLabel', br=FALSE)
            hwrite(queryArgs$jobID, p, br=TRUE)
        }
        
        hwrite("Raw gene-term sets text files: ", p,  br=FALSE)
        jobName <- gsub(.Platform$file.sep, "", folder)
        nRawMg <- nrow(feaResults$metagroups)
        hwrite(paste('[Mg', 1:nRawMg,']', sep=""), p, link=paste(folder,jobName,"_metagroup_", 1:nRawMg, ".txt", sep=""), table=FALSE, br=FALSE)
        hwrite('',p, br=TRUE)
    }
    
    hwrite("<br/>Formatted results from functional enrichment and clustering: ", p,  br=FALSE)
    curWd <- getwd()
    if(!grepl(paste("\\", .Platform$file.sep, "$",sep=""), curWd, fixed=FALSE)) curWd <- paste(curWd, .Platform$file.sep, sep="")
    fileName <- gsub(curWd, "", feaResults$fileName, fixed=TRUE)    
    hwrite(fileName,link=paste(folder, fileName, sep=""), p,  br=TRUE)
    

    ############################################
    ## Arguments
    hwrite("<br/>Arguments set for the query: ", p,  br=TRUE)        
    if(!is.null(queryArgs)) 
    {
        if("organism" %in% names(queryArgs))
        {
            data("organisms", envir = environment())
            organisms<- get("organisms", envir  = environment())
            hwrite("Organism: ", p, class='InfoLabel', br=FALSE) 
            hwrite(paste(queryArgs$organism, " (",organisms[queryArgs$organism,"Name"],")",sep=""), p, br=TRUE)
        }
        if("geneList" %in% names(queryArgs)) 
        {
            genesQuery <- strsplit(eval(queryArgs$geneList),", ")[[1]]
             hwrite("Genes: ", p, class='InfoLabel', br=FALSE)
             hwrite(paste("(", length(genesQuery), ") ", sep=""), p, br=FALSE)
             hwrite(paste(sort(genesQuery), sep="", collapse=", "), p, br=TRUE)
        }
    
        # Other arguments...    
        argumentsLeft <- names(queryArgs)[!names(queryArgs) %in% c("fun", "serverWeb", "geneList","organism", "argsWS", "...")] # argsWS se muestra muy raro
    
        # NOT eval: argsWS, eset, geneSets, genesUniverse, refPackage, geneID2GO
        for(arg in argumentsLeft)
        {
            hwrite(paste(arg, ": ",sep=""), p, class='InfoLabel', br=FALSE)
            hwrite(queryArgs[[arg]], p, br=TRUE)   
        }
    
        hwrite('</span>', p, br=TRUE) # Section end (contentQueryInfo)       
    }else{
        hwrite('Not available.', p, br=TRUE)
    }
}