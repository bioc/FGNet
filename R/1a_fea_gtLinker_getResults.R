### Functions: 
# getResults_gtLinker
# format_gtLinker (private)

# jobName <- "3907019_gtLinker_"
format_gtLinker <- function(jobName, organism)
{
    ############## READ Files ##############
    # Read metagroups
    globalMetagroups <- read.table(paste(jobName, "global_overview.txt", sep=""),skip=2, header=FALSE, sep="\t", quote = "")
    colnames(globalMetagroups) <- c("Size","Diameter","Similarity","Silhouette Width","Genes","nGenes","nref_list","pValue","Terms")
    nMetaGroups <- dim(globalMetagroups)[1]
    #write.table(globalMetagroups, file="...")
    
    # Read gene-term sets
    tablaGeneTermSets <- NULL
    for(mg in rownames(globalMetagroups)) 
    {
        metagrupo <-  read.table(paste(jobName, "metagroup_", mg, ".txt", sep=""),skip=1, header=FALSE, sep="\t", quote = "")
        tablaGeneTermSets <- rbind(tablaGeneTermSets,cbind(cbind(gtSet=paste(mg, sep=""), metagrupo))) #, cbind(gtSet=paste("Mg", mg,"_gtSet", 1:dim(metagrupo)[1], sep=""), metagrupo))
    }
    colnames(tablaGeneTermSets) <- c("Metagroup", "Genes","nGenes","nref_list","pValue","Terms")
    
    ############## Basic format (global metagroups and gsets) ###################
    # Remove brackets 
    for(col in c("nGenes", "nref_list"))
    {
        globalMetagroups[, col] <- as.numeric(sapply(strsplit(as.character(globalMetagroups[, "nGenes"]), split="(", fixed=TRUE), function(x) x[1]))
        tablaGeneTermSets[, col] <- as.numeric(sapply(strsplit(as.character(tablaGeneTermSets[, "nGenes"]), split="(", fixed=TRUE), function(x) x[1]))
    }
    
    # Format kegg, add org prefix...
    data("organisms", envir = environment())
    organisms <- get("organisms", envir  = environment())
    if(is.null(organism))
    {
        organism <- "map"    
    }else{
        organism <- organisms[capitalize(tolower(organism)), "keggPrefix"]
    }
    globalMetagroups$Terms <- gsub("Kegg:", paste("KEGG:", organism, sep=""), globalMetagroups$Terms)
    tablaGeneTermSets$Terms <- gsub("Kegg:", paste("KEGG:", organism, sep=""), tablaGeneTermSets$Terms)
    # Interpro
    globalMetagroups$Terms <- gsub("IPR", paste("IPR:", sep=""), globalMetagroups$Terms)
    tablaGeneTermSets$Terms <- gsub("IPR", paste("IPR:", sep=""), tablaGeneTermSets$Terms)
    
    ############## Format metagroups table ###################
    globalMetagroups <- cbind(Metagroup=1:nrow(globalMetagroups),globalMetagroups)    
    globalMetagroups$Terms <- sapply(globalMetagroups$Terms, function(gts)
            {
                paste(sapply(strsplit(gts,";", fixed=TRUE)[[1]], function(term)
                {
                    term <-strsplit(term, ":")[[1]]
                    capitalize(paste(paste(term[3:length(term)], collapse=" "), " (",paste(term[1:2], collapse=":"),")",sep=""))
                
                }),collapse=";")
            })
        
    ############## SAVE ###################
    fileNameMetagroups <- paste(jobName, "formatted_metagroups.txt", sep="")
    fileNameGeneTermSets <- paste(jobName, "formatted_gtSets.txt", sep="")
    
    if(!file.exists(file.path(fileNameMetagroups)))
    {
        write.table(globalMetagroups, file=fileNameMetagroups, quote=FALSE, row.names=FALSE, sep="\t")
        write.table(tablaGeneTermSets, file=fileNameGeneTermSets, quote=FALSE, row.names=FALSE, sep="\t")
        message(paste("Metagroups and Gene-term sets tables saved as ", fileNameMetagroups, " and ",fileNameGeneTermSets,".", sep=""))
    }
    
    ret <- list(metagroups=globalMetagroups, geneTermSets=tablaGeneTermSets, fileName=fileNameMetagroups)
}


# jobID <- 3907019
# results <- fea_gtLinker_getResults(jobID)
fea_gtLinker_getResults <- function(jobID=NULL, organism=NULL, jobName=NULL, alreadyDownloaded=FALSE, keepTrying=FALSE, serverWeb="http://gtlinker.cnb.csic.es", serverWS="http://gtlinker.cnb.csic.es:8182") 
{
    if(!loadInstPkg("RCurl")) stop("Package 'RCurl' is required to get the FEA results from GeneTerm Linker server.")

    # Check arguments
    if(is.numeric(jobID)) jobID <- as.character(jobID)
    if(!is.character(jobID) && !alreadyDownloaded) stop("jobID not valid.")
    if(is.null(jobID) && is.null(jobName)) stop("jobID or jobName are required.")
    
    if(!is.null(jobName) && !is.character(jobName)) stop("jobName should be character.")
    if(!is.logical(alreadyDownloaded)) stop("alreadyDownloaded should be TRUE or FALSE.")
    if(!is.logical(keepTrying)) stop("keepTrying should be TRUE or FALSE.")
    if(!is.character(serverWeb)) stop("Web server URL not valid.")
    if(!alreadyDownloaded && !url.exists(serverWeb)) stop("Either the web server URL is wrong or the server is down.")
    
    ###############################
    # Prepare jobName / folder
    if(is.null(jobName)) jobName <- paste(jobID,"_gtLinker", sep="")
    folder <- jobName
    jobName <- paste(jobName, "_", sep="")
    
    # Create folder
    if((!file.exists(file.path(folder))))
    {
        dir.create(file.path(folder))    
    }
    currWD <- getwd()
    setwd(folder)
    
    ##########################################################
    # Connect to server and download results
    ##########################################################
    if(!is.null(jobID) && !alreadyDownloaded)
    {
        jobReady <- FALSE
        firstAttempt <- TRUE
        count <- 0
        status_envelope_body <- paste('<status xmlns="urn:gtLinkerWS">',
                                       '<job_id xsi:type="xsd:string">',jobID ,'</job_id></status>', sep="") 
        
        ###############################
        # Check wether the job is ready and wait if necessary
        ###############################
        #     status(jobID)
        #       0: the job has finished succesfully
        #       1: the job has not finished yet (inc. the job does not exist!)
        #      -1: there has been an error 
        #     
        #     info(jobID) # Provides any error messages there might have been during the analisys.
        while(!jobReady)
        {
            reply <- SOAPQuery(status_envelope_body, serverWS)
                        
            # Ready
            if(reply$statusResponse == 0) 
            {
                jobReady <- TRUE
                if(!firstAttempt) Sys.sleep(5) # To allow the server producing the files 
            }
            
            # Error
            if(reply$statusResponse == -1)
            {
                info_envelope_body <- paste('<info xmlns="urn:gtLinkerWS">',
                                             '<job_id xsi:type="xsd:string">',jobID ,'</job_id></info>', sep="")
                reply <- SOAPQuery(info_envelope_body, serverWS)$infoResponse
                stop(paste("Message from server: ", reply, sep=""))
            }
            # if(reply == -1)  stop("Error connecting to the server.")
            
            # Still analyzing:
            if(!keepTrying) break()
            if(reply$statusResponse == 1) 
            {
                count <- count + 1
                firstAttempt <- FALSE
                if(count == 10) keepTrying <- FALSE
                Sys.sleep(10) 
            }
        }
        
        ###############################
        # Download results (jobID)
        ###############################
        if(!jobReady) 
        {
            warning("The analysis has not finished yet or the jobID does not exist.")
            return(NULL)
        } else  # jobReady
        {
            # Download files
            downAttempt <- tryCatch({
                # Download global_hash
                inputFileName <- paste(jobName,"global_overview.txt", sep="")  
                fileURL <- paste(serverWeb, "/jobs/",jobID,"/global_hash", sep="")
                if(!url.exists(fileURL)) warning("Server/URL does not seem to be correct, trying anyway...")
                download.file(url=fileURL, destfile=inputFileName)
                
                # Download Metagroups
                nMetaGroups <- nrow(read.table(inputFileName,skip=2, header=FALSE, sep="\t", quote = ""))
                for (mg in 1:nMetaGroups) download.file(quiet=TRUE, url=paste(serverWeb, "/jobs/",jobID,"/", mg-1,"_hash", sep=""), destfile=paste(jobName,"metagroup_", mg,".txt", sep=""))
                
                message(paste("Results downloaded to folder '", folder, "'", sep=""))
            }, error = function(e) {
                message(e)
                return(e$message)
            })
            if(!is.null(downAttempt)) return(downAttempt)
        }
    }# End download
    
    queryArgs <- list(queryArgsAsCharacter(match.call()))
    ret <- c(queryArgs=queryArgs, format_gtLinker(jobName, organism))
    setwd(currWD)
    
    invisible(ret)
}
