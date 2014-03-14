# jobID <- 3907019
# results <- getResults_gtLinker(jobID)

getResults_gtLinker <- function(jobID=NULL,  path=getwd(), jobName="", alreadyDownloaded=FALSE, keepTrying=FALSE, serverWeb="http://gtlinker.cnb.csic.es", serverWS="http://gtlinker.cnb.csic.es:8182")
{
	# Libraries
	#	if (!library(igraph, logical.return=TRUE)) stop("Library igraph is required to plot the networks.")
	
	# Check arguments
	if(is.numeric(jobID)) jobID <- as.character(jobID)
	if(!is.character(jobID) && !alreadyDownloaded) stop("jobID not valid.")
	
	if(!file.exists(path)) stop("The given path does not exist.")
	if(!is.character(jobName)) stop("jobName should be character.")
	if(!is.logical(alreadyDownloaded)) stop("alreadyDownloaded should be TRUE or FALSE.")
	if(!is.logical(keepTrying)) stop("keepTrying should be TRUE or FALSE.")
	if(!is.character(serverWeb)) stop("Web server URL not valid.")
	if(!url.exists(serverWeb)) stop("Either the web server URL is wrong or the server is down.")
	
	# Create folder
	if(substring(path, first=nchar(path), last=nchar(path)) != "/") path <- paste(path, "/", sep="")
	folder <- path
	if(jobName != "") 
	{
		folder <- paste(path, jobName, "/", sep="")	
		jobName <- paste(jobName, "_", sep="")
	}
	
	if ((!file.exists(folder)) && (!folder==""))
	{
		dir.create(file.path(folder))		
	}
	
	###############################
	# Check wether the job is ready
	###############################
	# 	status(jobID)
	# 	  0: the job has finished succesfully
	# 	  1: the job has not finished yet (inc. the job does not exist!)
	# 	 -1: there has been an error 
	# 	
	# 	info(jobID) # Provides any error messages there might have been during the analisys.
	
	jobReady <- FALSE
	firstAttempt <- TRUE
	count <- 0
	if(!is.null(jobID))
	{
		status_envelope_body <-	 paste('<status xmlns="urn:gtLinkerWS">',
																	 '<job_id xsi:type="xsd:string">',jobID ,'</job_id></status>', sep="") 
		
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
				info_envelope_body <-	 paste('<info xmlns="urn:gtLinkerWS">',
																		 '<job_id xsi:type="xsd:string">',jobID ,'</job_id></info>', sep="")
				reply <- SOAPQuery(info_envelope_body, serverWS)$infoResponse
				stop(paste("Message from server: ", reply, sep=""))
			}
			
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
	}								
	
	###############################
	# Download results (jobID)
	###############################
	if(!jobReady) 
	{
		warning("The analysis has not finished yet or the jobID does not exist.")
		return(NULL)
	}
	if(jobReady)
	{
		downAttempt <- tryCatch({
			# Global
			inputFileName <- paste(folder, jobName,"global_overview.txt", sep="")
			if(!alreadyDownloaded) 
			{
				fileURL <- paste(serverWeb, "/jobs/",jobID,"/global_hash", sep="")
				if(!url.exists(fileURL)) warning("Server/URL does not seem to be correct, trying anyway...")
				download.file(url=fileURL, destfile=inputFileName)
			}
			globalMetagroups <- read.table( paste(folder, jobName, "global_overview.txt", sep="") ,skip=2, header=F, sep="\t", quote = "")
			colnames(globalMetagroups) <- c("Size","Diameter","Similarity","Silhouette Width","Genes","nGenes","nref_list","pValue","Terms")
			nMetaGroups <- dim(globalMetagroups)[1]
			
			globalMetagroups[, "nGenes"] <- as.numeric(sapply(strsplit(as.character(globalMetagroups[, "nGenes"]), split="(", fixed=TRUE), function(x) x[1]))
			globalMetagroups[, "nref_list"] <- sapply(strsplit(as.character(globalMetagroups[, "nref_list"]), split="(", fixed=TRUE), function(x) x[1])
			
			# Download Metagroups
			if(!alreadyDownloaded) for (mg in 1:nMetaGroups) download.file(quiet=TRUE, url=paste(serverWeb, "/jobs/",jobID,"/", mg-1,"_hash", sep=""), destfile=paste(folder, jobName,"metagroup_", mg,".txt", sep=""))
			
			# Leer termSets
			tablaGeneTermSets <- NULL
			for(mg in rownames(globalMetagroups)) 
			{
				metagrupo <-  read.table(paste(folder, jobName, "metagroup_", mg, ".txt", sep=""),skip=1, header=F, sep="\t", quote = "")
				tablaGeneTermSets <- rbind(tablaGeneTermSets,cbind(cbind(gtSet=paste(mg, sep=""), metagrupo))) #, cbind(gtSet=paste("Mg", mg,"_gtSet", 1:dim(metagrupo)[1], sep=""), metagrupo))
			}
			colnames(tablaGeneTermSets) <- c("Metagroup", "Genes","nGenes","nref_list","pValue","Terms")
			
			if(!alreadyDownloaded) message(paste("Results downloaded to ", folder, sep=""))
			return(list(metagroups=globalMetagroups, geneTermSets=tablaGeneTermSets, fileName=inputFileName))
		}, error = function(e) {
			FALSE
		})
		if(downAttempt == FALSE) 	
		{
			if(!jobReady) message("The analysis has not finished yet or the jobID does not exist.")
		} else {
			return (downAttempt)
		}
	}
}
