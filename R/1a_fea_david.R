
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

# geneIdType <- "GENE_SYMBOL"
# annotations <- c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO")
# geneList <- c("A2M", "ABL1", "APBA1", "APBB1", "APLP1", "APLP2", "APOE", "APP", "ATOH1", "BRCA1", "BRCA2", "CDK5R1", "CDK5", "CDK5R2", "DAB1", "DLL1", "DNMT1", "EGFR", "ERBB2", "ETS1", "FOS", "FYN", "GLI1", "GLI2", "GLI3", "JAG1", "KIT", "LRP1", "LRP8", "MAPT", "MYC", "NOTCH1", "NRAS", "PAX2", "PAX3", "PSEN1", "PSEN2", "PTCH1", "RELN", "ROBO1", "SHC1", "SHH", "SMO", "SRC", "TGFB1", "TP53", "VLDLR", "WNT1", "WNT2", "WNT3")

# API: david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
# Maximum gene ids: 4000
fea_david <- function(geneList, geneIdType="ENSEMBL_GENE_ID", geneLabels=NULL, annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), email=NULL, argsWS=c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L), jobName=NULL, downloadFile=TRUE)
{    
    # Check arguments
    if(!is.character(geneList)) stop("geneList should be a character vector.")
    
    # argsWS
    if(is.null(email))
    {
        argsWS <- NULL
    }else{
        argsNames <- c("overlap", "initialSeed", "finalSeed", "linkage", "kappa")
        if(!is.numeric(argsWS)) stop("argsWS should contain the arguments to pass to the web service.")
        if(any(!names(argsWS) %in% argsNames)) stop(paste("Web service arguments should be: ", paste(argsNames, collapse=", "), ".", sep=""))
        if(any(!argsNames %in% names(argsWS)))
        {
            defaultArgs <- c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L)
            argsNames <- c(argsWS, defaultArgs[argsNames[which(!argsNames %in% names(argsWS))]])
        }  
    }

    # Format arguments:
    queryArgs <- list(queryArgsAsCharacter(match.call()))
    
    # Query:
    if(!is.null(email)) # WebService
    {        
        if(!loadInstPkg("RDAVIDWebService")) stop("Package RDAVIDWebService is required to query DAVID through the webserver. Install the package or set email=NULL to query DAVID through the web API.")
        
        randomNumber <- sample(100000:999999, 1)
        if(is.null(jobName)) jobName <- paste(randomNumber, "_DAVID", sep="")
        downloadFileName <- paste(tempdir(), .Platform$file.sep, jobName, ".txt", sep="")
        
        # Connect to DAVID        
        tryCatch( 
        {
            davidConnection <- DAVIDWebService$new(email=email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
        }, warning = function (w)
        {
            if(grepl("SSL", w)) w <- paste("If you are seeing this error you might need to install a certificate to connect to DAVID web service.\n See the instructions in the forum: \nhttps://support.bioconductor.org/p/70090/#72226", w, sep="\n")
            errorMsgDavid(w)
        })
                
        # Upload gene list
        result <- tryCatch( 
            {
                addList(davidConnection, geneList, idType=geneIdType, listName=paste("List_", randomNumber, sep=""), listType="Gene")
            }, error = function(e) 
            {
                if(grep("idType", e$message)) 
                {
                    idsFile <- "DAVID_GeneIdTypes.txt"
                    write.table(getIdTypes(davidConnection), file=idsFile, row.names=FALSE, col.names=FALSE, quote=FALSE)
                    
                    #errorMsg <- simpleError(message=paste("Gene ID type not valid. Available IDs: ", paste(getIdTypes(davidConnection), collapse=", "), ".",sep=""), call=functionCall)
                    errorMsg <- simpleError(message=paste("Gene ID type not valid. Available ID types were saved in the file ", idsFile, sep=""), call=match.call())
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
                
                errorMsg <- simpleError(message=paste("Some of the annotation categories are not valid. Available ones were saved in the file ", annotFile, sep=""))
                stop(errorMsg)
            } else {
                errorMsgDavid(e)
            }
        }
        )
        
        # Request & save clustering report
        getClusterReportFile(davidConnection, type="Term", fileName=downloadFileName, 
                                                 overlap=argsWS["overlap"], initialSeed=argsWS["initialSeed"], finalSeed=argsWS["finalSeed"], linkage=argsWS["linkage"], kappa=argsWS["kappa"])
        moveFile <- TRUE
        
    }else
    {
        ## Query through web API
        if(!loadInstPkg("RCurl")) stop("We highly recommend to provide a registered email to query DAVID through the Web Service.\nOtherwise, please install package 'RCurl' to continue.")
        
        if(length(geneList)>400) stop("Maximum number of genes: 400. For more genes register at DAVID Web Service (see help for details).")
        geneList <- paste(geneList, collapse=",")
        annotations <- paste(annotations, collapse=",")
        
        tool <- "term2term"
        queryUrl <- paste("https://david.ncifcrf.gov/api.jsp?type=", geneIdType, "&ids=", geneList, "&tool=", tool, "&annot=", annotations, sep="") # Do not change order
        
        # URL has a charactor size limitation (<= 2048 characters in total), i.e., the very large gene list may not be able to completely passed by URL.
        if(nchar(queryUrl)>2048) stop("Query url too long.")
                
        curlHandle <- getCurlHandle(cookiefile = "CurlHandleCookie.txt")        
        reply <- getURL(queryUrl, curl = curlHandle, ssl.verifypeer = FALSE)

        replyRowids <- getContent(reply, attribute = 'document.apiForm.rowids.value="')
        replyAnnot <- getContent(reply, attribute = 'document.apiForm.annot.value="')
        
        getURL <- paste("https://david.ncifcrf.gov/term2term.jsp?rowids=", replyRowids, "&annot=", replyAnnot, sep="")
        if(nchar(getURL) < 2048)
        {
            finalReply <- getURL(getURL, curl = curlHandle, ssl.verifypeer = FALSE)
        }else
        {    
            finalReply <- tryCatch({
                postForm("https://david.ncifcrf.gov/term2term.jsp",
                                 curl = curlHandle,
                                 rowids = replyRowids,
                                 annot = replyAnnot)
            }, error = function(e) {
                FALSE
            })
            
            if(finalReply != FALSE) warning("Query URL too long, default annotations might be used.")  
            if(finalReply == FALSE) stop("Query URL too long, try with less genes or use the Web Server (provide email).")  
        }
    
        downloadFileName <- NULL
        if(finalReply != FALSE)
        {
            downloadFileName <- getContent(finalReply, attribute = '<a href="data/download/')[1]
        }
        
        if(grepl(".txt", downloadFileName))
        {
            downloadFileName <- paste("https://david.ncifcrf.gov/data/download/", downloadFileName, sep="")
            moveFile <- FALSE
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
print(downloadFileName)
print(jobName)
print(moveFile)

    ret <- format_david(downloadFileName, jobName=jobName, geneLabels=geneLabels, moveFile=moveFile, downloadFile=downloadFile)
    invisible(c(queryArgs=queryArgs,ret))
}

