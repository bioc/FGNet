# Workflow wrapper, generates a Report in HTML.

report_gtLinker <- function(genes=NULL, organism="Hs", annotations=c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs"), minSupport=4, jobID=NULL, alreadyDownloaded=FALSE, path=getwd(), jobName=NULL, threshold=0, serverWeb="http://gtlinker.cnb.csic.es", serverWS="http://gtlinker.cnb.csic.es:8182")
{
	return( fnReport(tool="gtLinker", organism=organism, genes=genes, annotations=annotations, minSupport=minSupport, jobID=jobID, alreadyDownloaded=alreadyDownloaded, path=path, jobName=jobName, threshold=threshold, serverWeb=serverWeb, serverWS=serverWS))
}
report_david <- function(genes=NULL, geneIdType="ENSEMBL_GENE_ID", annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), email=NULL,  argsWS = c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L), inputFileLocation=NULL, path=getwd(), jobName=NULL, threshold=0, geneLabels=NULL)
{
	return( fnReport(tool="David", genes=genes, geneIdType=geneIdType, annotations=annotations, email=email, inputFileLocation=inputFileLocation, alreadyDownloaded=FALSE, path=path, jobName=jobName, threshold=threshold, argsWS = argsWS, geneLabels=geneLabels))
}

### Common arguments
# tool
# genes
# annotations (Different annotation name)
# alreadyDownloaded
# path
# jobName
# threshold
# 
### Arguments only for gtLinker:
# organism
# minSupport
# jobID
# serverWeb
# serverWS
# 
### Arguments only for David:
# geneIdType
# inputFileLocation

fnReport <- function (tool, organism=NULL, geneIdType=NULL, genes=NULL, annotations=NULL, email=NULL, minSupport=3, jobID=NULL, inputFileLocation=NULL, alreadyDownloaded=FALSE, path=NULL, jobName=NULL,  threshold=0, serverWeb=NULL, serverWS=NULL, argsWS = NULL, geneLabels=geneLabels)
{
	#####################################################################################################
	####################################   Check arguments   ############################################
	if(tool == "gtLinker")
	{
		if(is.null(jobID))
		{
			if(is.null(genes) && !alreadyDownloaded) stop("If the analisys is not already downloaded, either a jobID or a query (organism, genes and annotations) should be provided.")
						
			if(!is.character(organism)) stop("Organism should be either 'Hs' (Homo sapiens) or 'Sc' (Saccharomyces cerevisiae). ")
				organism <- tolower(organism) 
				if (!(organism %in% c("sc","hs"))) stop("Organism should be either 'Hs' (Homo sapiens) or 'Sc' (Saccharomyces cerevisiae).")

			if(!is.numeric(minSupport)) stop("minSupport should be a number.")		
				minSupport <- as.integer(minSupport)
				if(minSupport<2 || minSupport>100) stop("minSupport should be more than 1.")		
		} else 
		{
			if(!is.null(genes)) stop("Either a jobID *OR* a query (organism, genes and annotations) should be provided.")
			if(is.numeric(jobID)) jobID <- as.character(jobID)
			if(!is.character(jobID)) stop("jobID not valid.")
		}
		
		if(!is.character(serverWS)) stop("Webservice server url not valid.")
		if(!is.character(serverWeb)) stop("Web server url not valid.")
	} else 
	{
		if(tool == "David") 
		{
			if(!is.character(geneIdType)) stop("geneIdType not valid.") # List of available ones in David API
			if(!is.null(inputFileLocation) && !is.character(inputFileLocation)) stop("inputFileLocation not valid.")
			if(is.null(inputFileLocation) && is.null(genes))  stop("Either a inputFileLocation a query (organism, genes and annotations) should be provided.")
			if(!is.null(inputFileLocation) && !is.null(genes))  stop("Either an inputFileLocation *OR* a query (genes and annotations) should be provided.")
		}else
		{
			stop("'tool' not valid.")
		}
	}
	
	if(!is.null(genes) && !is.character(genes)) stop("Genes should be a character vector containing the names of the genes.")
	
	if(tool == "gtLinker")
	{
			allowedAnnotations <- c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs")
			if (any(!tolower(annotations) %in% tolower(allowedAnnotations))) warning(paste("Some of the provided annotations are not recognized. Example of available annotations: ", paste(allowedAnnotations, collapse=", "), sep=""))	
	}
	
	if(!is.character(annotations)) stop("Annotations should be a charater vector.")

	if(!is.logical(alreadyDownloaded)) stop("alreadyDownloaded should be logical.")
	if(!file.exists(path)) stop("The given path does not exist.")
	if(!is.null(jobName) && !is.character(jobName)) stop("jobName should be character.")
	if(!is.numeric(threshold)) stop("threshold should be numeric")
	
	
	if(!(is.null(serverWeb) && is.null(serverWS)) && !is.null(genes))
	{
		if((is.null(serverWeb) + is.null(serverWS)) == 1) 
		{
			warning("'serverWS' and 'serverWeb' doesn't seem to match.")
			
		}else {
			ws <- strsplit(sub("http://","", serverWS), split=":")[[1]][1]
			web <- strsplit(sub("http://","", serverWeb), split=":")[[1]][1]
			if(! ws %in% web) warning("The servers for the webservice and the web doesn't seem to match.")
		}
	}
	
	#####################################################################################################
	################################## Query and get Results ############################################
	# gtLinker
	if(tool == "gtLinker")
	{
		grType <- "Metagroup"
		metagroupAttributeName <- "Silhouette Width"
		
		if(alreadyDownloaded && (is.null(jobID) && is.null(jobName))) stop("If the analysis is already downloaded, please provide the jobID or jobName (folder).") #Enough??
				
		if(is.null(jobID) && !alreadyDownloaded)
		{
			jobID <- query_gtLinker(genes=genes, organism=organism, annotations=annotations, minSupport=minSupport, serverWS=serverWS)
			if (jobID == "-1") stop ("There was an error during the query to GeneTerm Linker. Please, try again.")
		}else{
			annotations <- NULL
			organism  <- NULL
			genes <- NULL
		}
		if(is.null(jobName)) jobName <- paste("FunctionalNetwork_", jobID, "_thr", threshold, sep="")
		results <- getResults_gtLinker(jobID=jobID, keepTrying=TRUE, path=path, jobName=jobName, alreadyDownloaded=alreadyDownloaded, serverWeb=serverWeb)		
		if(is.null(results)) stop("GtLinked does not return any results.")
	}
	
	# David
	if(tool == "David")
	{
		metagroupAttributeName <- "EnrichmentScore"
		serverWeb <- "http://david.abcc.ncifcrf.gov/summary.jsp"
	
		if(is.null(inputFileLocation))
		{
			inputFileLocation <- query_david(genes=genes, geneIdType=geneIdType, annotations=annotations, email=email, argsWS=argsWS)
			fileProvided <- FALSE
		}else{
			fileProvided <- TRUE
		}
		if(is.null(jobName))
		{
			if(grepl("http://", inputFileLocation))
			{
				prefix <- "FunctionalNetwork"
				sufix <- strsplit(inputFileLocation, "/download/", fixed="TRUE")[[1]][2]
				sufix <- paste(substr(sufix, start=5, stop=7), substr(sufix, start=nchar(sufix)-7, stop=nchar(sufix)-4), sep="")
			}else # Local file
			{
				annotations <- NULL
				first <- gregexpr("./.", inputFileLocation)[[1]]
				first <- first[length(first)]+2
			
				last <- gregexpr(".txt", inputFileLocation)[[1]][1]-1
				prefix <- substring(inputFileLocation, first=first, last=last) # Location (Folder)
			}
			sufix <- paste("_thr", threshold, sep="")
			if(!fileProvided)
			{
				jobName <- paste(prefix, sufix, sep="")
			}else{
				jobName <- ""
			}
		}
		if(!is.null(geneLabels) && is.null(names(geneLabels)))
		{
			warning("geneLabels not valid.")
			geneLabels <- NULL
		}
		results <- getResults_david(inputFileLocation=inputFileLocation, path=path, jobName=jobName, geneLabels=geneLabels)
		if(is.null(annotations)) annotations <- levels(results$geneTermSets$Category)
	}
	
	# Common to both tools
	attribute <- results[[1]][,metagroupAttributeName, drop=FALSE]   # $metagroups or $cluster (first in list)
	tables <- adjMatrix(results$geneTermSets, attribute = attribute, threshold=threshold)

	#####################################################################################################
	####################################   Generate  HTML    ############################################
	
	htmlFileName <- paste(path, "/", jobName, ".html", sep="")
	
	createHtml(htmlFileName=htmlFileName, results=results, tables=tables, metagroupAttributeName=metagroupAttributeName, threshold=threshold,
						 organism=organism, annotations=annotations, genes=genes, serverWeb=serverWeb, jobID=jobID)

	#return (jobID)
}

