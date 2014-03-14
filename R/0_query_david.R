
# Reads content until "
getContent <- function(reply, attribute)
{
	start <- gregexpr(attribute, reply)[[1]]
	start <- start[1] + attr(start, "match.length")
	end <- regexpr('\"', substring(reply, start, nchar(reply))) + start -2
	
	ret <- substr(reply, start=start, stop=end)
	return(ret)
}

errorMsgDavid <- function(errorMsg)
{
	stop(paste("Message from DAVID server: ", errorMsg$message, sep=""), call.=FALSE)
}

#geneIdType <- "GENE_SYMBOL"
#annotations <- c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO")
#genes <- c("A2M", "ABL1", "APBA1", "APBB1", "APLP1", "APLP2", "APOE", "APP", "ATOH1", "BRCA1", "BRCA2", "CDK5R1", "CDK5", "CDK5R2", "DAB1", "DLL1", "DNMT1", "EGFR", "ERBB2", "ETS1", "FOS", "FYN", "GLI1", "GLI2", "GLI3", "JAG1", "KIT", "LRP1", "LRP8", "MAPT", "MYC", "NOTCH1", "NRAS", "PAX2", "PAX3", "PSEN1", "PSEN2", "PTCH1", "RELN", "ROBO1", "SHC1", "SHH", "SMO", "SRC", "TGFB1", "TP53", "VLDLR", "WNT1", "WNT2", "WNT3")
#"http://david.abcc.ncifcrf.gov/api.jsp?type=GENE_SYMBOL&ids=SRC,TGFB1,TP53,VLDLR,WNT1,WNT2,WNT3&tool=term2term&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"
#fileName<-query_david(genes=c("A2M", "ABL1", "APBA1", "APBB1", "APLP1", "APLP2", "APOE", "APP", "ATOH1", "BRCA1", "BRCA2", "CDK5R1", "CDK5", "CDK5R2", "DAB1", "DLL1", "DNMT1", "EGFR", "ERBB2", "ETS1", "FOS", "FYN", "GLI1", "GLI2", "GLI3", "JAG1", "KIT", "LRP1", "LRP8", "MAPT", "MYC", "NOTCH1", "NRAS", "PAX2", "PAX3", "PSEN1", "PSEN2", "PTCH1", "RELN", "ROBO1", "SHC1", "SHH", "SMO", "SRC", "TGFB1", "TP53", "VLDLR", "WNT1", "WNT2", "WNT3"), annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), geneIdType="GENE_SYMBOL")

# API: david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
# Maximum gene ids: 400
query_david <- function(genes, geneIdType="ENSEMBL_GENE_ID", annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), email=NULL, argsWS = c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L))
{	
	# Check arguments
	if(!is.character(genes)) stop("genes should be a character vector.")
	
	functionCall <- match.call()
	if(!is.null(email))
	{
		if(!library(RDAVIDWebService, logical.return=TRUE, quietly=TRUE)) stop("Package RDAVIDWebService is required t query DAVID through the webserver. Install the package or set email=NULL to query DAVID through the web API.")
		randomNumber <- sample(100:1000, 1)
		downloadFile <- paste("DavidClustering_", randomNumber, ".txt", sep="")
		
		# Check argsWS
		#argsWS = list(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L)
		argsNames <- c("overlap", "initialSeed", "finalSeed", "linkage", "kappa")
		if(!is.numeric(argsWS)) stop("argsWS should contain the arguments to pass to the web service.")
		if(any(!names(argsWS) %in% argsNames)) stop(paste("Web service arguments should be: ", paste(argsNames, collapse=", "), ".", sep=""))
		if(any(!argsNames %in% names(argsWS)))
		{
			defaultArgs <- c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L)
			argsNames <- c(argsWS, defaultArgs[argsNames[which(!argsNames %in% names(argsWS))]])
		}
		
		# Connect to DAVID		
		tryCatch( 
{
	davidConnection <- DAVIDWebService$new(email=email)
}, warning = function (w)
{
	errorMsgDavid(w)
})
		
		# Upload gene list
		tryCatch( 
{
	result <- addList(davidConnection, genes, idType=geneIdType, listName=paste("List_", randomNumber, sep=""), listType="Gene")
}, error = function(e) 
{
	if(grep("idType", e$message)) 
	{
		idsFile <- "DAVID_GeneIdTypes.txt"
		write.table(getIdTypes(davidConnection), file=idsFile, row.names=FALSE, col.names=FALSE, quote=FALSE)
		
		#errorMsg <- simpleError(message=paste("Gene ID type not valid. Available IDs: ", paste(getIdTypes(davidConnection), collapse=", "), ".",sep=""), call=functionCall)
		errorMsg <- simpleError(message=paste("Gene ID type not valid. Available ID types were saved in the file ", idsFile, sep=""), call=functionCall)
		stop(errorMsg)
	} else {
		errorMsgDavid(e)
	}
}
		)
		if(length(result$unmappedIds) > 0) message(paste("Unmapped IDs: ", paste(result$unmappedIds, collapse=", "), ".", sep=""))	
		
		# Set annotations
		tryCatch( 
{
	setAnnotationCategories(davidConnection, categories=annotations)					
}, error = function(e) 
{
	if(grep("categories", e$message)) 
	{
		annotFile <- "DAVID_AnnotationCategoryNames.txt"
		write.table(getAllAnnotationCategoryNames(davidConnection), file=annotFile, row.names=FALSE, col.names=FALSE, quote=FALSE)
		
		errorMsg <- simpleError(message=paste("Some of the annotation categories are not valid. Available ones were saved in the file ", annotFile, sep=""), call=functionCall)
		stop(errorMsg)
	} else {
		errorMsgDavid(e)
	}
}
		)
		
		# Check species
		# specie <- getSpecieNames(davidConnection)
		# if(length(specie)>1)
		# {
		# 	if(org=="Hs") specie <- which(specie == "Homo sapiens(49)")
		# 	if(org=="Sc") specie <- which(specie == "Saccharomyces cerevisiae(58)")
		# 	setCurrentSpecies(davidConnection, species=specie)
		# }
		
		# Request & save clustering report
		getClusterReportFile(davidConnection, type="Term", fileName=downloadFile,
												 overlap=argsWS["overlap"], initialSeed=argsWS["initialSeed"], finalSeed=argsWS["finalSeed"], linkage=argsWS["linkage"], kappa=argsWS["kappa"])
		
	}else
	{
		## Query through web API
		
		if(length(genes)>400) stop("Maximum number of genes: 400")
		genes <- paste(genes, collapse=",")
		annotations <- paste(annotations, collapse=",")
		
		tool <- "term2term"
		queryUrl <- paste("http://david.abcc.ncifcrf.gov/api.jsp?type=", geneIdType, "&ids=", genes, "&tool=", tool, "&annot=", annotations, sep="") # Do not change order
		# URL has a charactor size limitation (<= 2048 characters in total), i.e., the very large gene list may not be able to completely passed by URL.
		if(nchar(queryUrl)>2048) stop("Query url too long.")
		
		writeChar(queryUrl, "queryUrl.txt")
		
		curlHandle <- getCurlHandle(cookiefile = "CurlHandleCookie.txt")
		reply <- getURL(queryUrl, curl = curlHandle)
		#writeChar(reply, "reply.txt")
		
		replyRowids <- getContent(reply, attribute = 'document.apiForm.rowids.value="')
		replyAnnot <- getContent(reply, attribute = 'document.apiForm.annot.value="')
		
		getURL <- paste("http://david.abcc.ncifcrf.gov/term2term.jsp?rowids=", replyRowids, "&annot=", replyAnnot, sep="")
		if(nchar(getURL) < 2048)
		{
			finalReply <- getURL(getURL, curl = curlHandle)
		}else
		{    
			finalReply <- tryCatch({
				postForm("http://david.abcc.ncifcrf.gov/term2term.jsp",
								 curl = curlHandle,
								 rowids = replyRowids,
								 annot = replyAnnot)
			}, error = function(e) {
				FALSE
			})
			
			if(finalReply != FALSE) warning("Query URL too long, default annotations might be used.")  
			if(finalReply == FALSE) stop("Query URL too long, try with less genes or use the Web Server (provide email).")  
		}
		
		downloadFile <- NULL
		if(finalReply != FALSE)
		{
			downloadFile <- getContent(finalReply, attribute = '<a href="data/download/')[1]
		}
		
		
		if(grepl(".txt", downloadFile))
		{
			downloadFile <- paste("http://david.abcc.ncifcrf.gov/data/download/", downloadFile, sep="")
		}else 
		{
			errorMsg <- getContent(finalReply, attribute = '<div class="error">')[1]
			errorMsg <- substr(errorMsg, start=1, stop=regexec("</div>", errorMsg)[[1]][1]-1)
			if(errorMsg == "Warning: Your list contains more than 3000 genes.  Please select a list with less than 3000 genes to use this tool.")
			{
				stop("Too many genes for David API. Try with diferent gene ID or executing it manualy through the web and downloading the results as txt.")
			} 
			if(errorMsg == "")
			{
				stop("Error in query to David. Please, check if you are providing the right ID type. ")
			}
			stop(paste("Message from David: ", errorMsg, sep=""))
		}
	} 
	return (downloadFile)
}
