
### Index:
# actDeactGenes
# processGenesField
# processExpressionField
# submitQuery
# loadFileDialog
# readFEAresults

#argsList: geneListFrame, esetBox
actDeactGenes <- function(pointer1, pointer2, newPage, argsList) 
{
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    # tabsTool$GetCurrentPage()  (OLD page)
    curPage <- rownames(FEA_tools)[which(FEA_tools[,"ID"]==newPage)]
    if(curPage %in% c("GeneTerm Linker", "topGO")) #"DAVID", 
    {
        RGtk2::gtkWidgetSetSensitive(argsList$geneListFrame, TRUE)   
        argsList$geneListFrame$visible <- TRUE
        
        RGtk2::gtkWidgetSetSensitive(argsList$esetBox, FALSE)
        argsList$esetBox$visible <- FALSE
    }
    if(curPage %in% c("gage"))
    {
        RGtk2::gtkWidgetSetSensitive(argsList$geneListFrame, FALSE)   
        argsList$geneListFrame$visible <- FALSE
        
        RGtk2::gtkWidgetSetSensitive(argsList$esetBox, TRUE)   
        argsList$esetBox$visible <- TRUE
    }
    if(curPage %in% c("Imported text file"))
    {
        RGtk2::gtkWidgetSetSensitive(argsList$geneListFrame, FALSE)   
        RGtk2::gtkWidgetSetSensitive(argsList$esetBox, FALSE)   
    }
}

processGenesField <- function(genesText)
{
    sepChar <- "\n"
    genesBuffer <- RGtk2::gtkTextViewGetBuffer(genesText)
    
    geneList <- RGtk2::gtkTextBufferGetText(genesBuffer, RGtk2::gtkTextBufferGetStartIter(genesBuffer)$iter, RGtk2::gtkTextBufferGetEndIter(genesBuffer)$iter, include.hidden.chars = TRUE)
    geneList <- unlist(strsplit(geneList, sepChar))
    
    # Clear whitespace
    # geneList <- unlist(strsplit(geneList, " "))
    geneList <- sub(" ", "", geneList)
    
    # Clear empty slots
    geneList <- geneList[geneList!=""]
    
    return(geneList)
}

processExpressionField <- function(expressionText)
{
    genesBuffer <- RGtk2::gtkTextViewGetBuffer(expressionText)
    geneList <- RGtk2::gtkTextBufferGetText(genesBuffer, RGtk2::gtkTextBufferGetStartIter(genesBuffer)$iter, RGtk2::gtkTextBufferGetEndIter(genesBuffer)$iter, include.hidden.chars = TRUE)
    geneList <- unlist(strsplit(geneList, "\n"))
    
    # Clear whitespace
    geneList <- sub(" ", "", geneList)
    
    # Clear empty slots
    geneList <- geneList[geneList!=""]
    
    # Get geneExpr
    geneList <- do.call(rbind,sapply(geneList, function(x) { 
                                                x <- gsub(" ", "\t",x)
                                                strsplit(x, "\t")
                                                }))
    geneList <- setNames(as.numeric(geneList[,2]), geneList[,1])
    
    if(length(geneList)==0) geneList <- NULL
    
    return(geneList)
}

############################################################################################################
# Used by: GUI_tabNetwork_common, GUI_tabFEA_other
# Text field is selected based on button name
# argsList= parentWindow +
# buttonLoadGenes: genesText
# buttonSelectExprFile: expressionText
# buttonSelectFeaResults: feaResultsText, comboFeaTool, expressionText, alreadyDownloadedCheck
# buttonImportFeaResults: feaResultsImportText
# esetButton: esetTxt
# geneSetsButton: geneSetsTxt$setText

loadFileDialog <- function(button, argsList) # button: Not used, but required for signalConnect.
{ 
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    fileName <- ""

    ### File dialog:
    dialog <- RGtk2::gtkFileChooserDialog("Choose File", argsList$parentWindow, "open",
                                   "gtk-open", RGtk2::GtkResponseType["accept"],
                                   "gtk-cancel", RGtk2::GtkResponseType["cancel"])
    
    if(dialog$run() == RGtk2::GtkResponseType["accept"]) 
    {
        fileName <- dialog$getFilename()
        
        # Replace backlash: fileName <- "C:\\blablabla\\blablabla"
        if(grepl("^[a-zA-Z]{1}:\\\\", fileName, fixed=FALSE)) fileName <- gsub("\\", .Platform$file.sep, fileName, fixed=TRUE)
        
        if(button$name == "buttonLoadGenes") 
        {
            ### Read genes...
            geneList <- as.character(unlist(read.table(fileName, sep="\t")[,1]))
            geneList <- paste(geneList,collapse="\n")
            # TO DO: Add file options?
            
            # Add to genes field
            genesBuffer <- RGtk2::gtkTextBufferNew()
            RGtk2::gtkTextBufferSetText(genesBuffer, geneList)
            RGtk2::gtkTextViewSetBuffer(argsList$genesText, genesBuffer)
        }
        if(button$name == "buttonSelectExprFile") 
        {
            ### Read genes...
            geneExpr <- read.table(fileName, sep="\t")
            geneExpr <- paste(paste(geneExpr[,1], geneExpr[,2], sep="\t"), collapse="\n")
            # TO DO: Add file options?
            
            # Add to genes field
            genesBuffer <- RGtk2::gtkTextBufferNew()
            RGtk2::gtkTextBufferSetText(genesBuffer, geneExpr)
            RGtk2::gtkTextViewSetBuffer(argsList$expressionText, genesBuffer)
        }
    
        if(button$name == "buttonSelectFeaResults") 
        {
            argsList$feaResultsText$setText(fileName)
    
            # For GeneTerm Linker
            if(grepl("global_overview.txt",fileName) || grepl("metagroup",fileName))
            {
                data("FEA_tools", envir = environment())
                FEA_tools<- get("FEA_tools", envir  = environment())
        
                argsList$alreadyDownloadedCheck$active <- TRUE
                argsList$comboFeaTool$setActive(as.numeric(FEA_tools["GeneTerm Linker", "ID"]))
            }
            
            ## Cambiar a changed feaResultsText? cuidado ftlinker (jobID)
            feaResults <- readFEAresults(feaResultsText=argsList$feaResultsText, comboFeaTool=argsList$comboFeaTool, expressionText=argsList$expressionText,
                                                                  serverWebReportText=argsList$serverWebReportText, alreadyDownloadedCheck=argsList$alreadyDownloadedCheck, parentWindow=argsList$parentWindow)
        }
        
        if(button$name == "buttonImportFeaResults")
        {
            argsList$fields$feaResultsImportText$setText(fileName)
            processFile(fileName, fieldsList=argsList$fields)
        }
    
        if(button$name == "esetButton")
        {
            argsList$esetTxt$setText(fileName)
            loadEset(fileName, refSamplesBox=argsList$refSamplesBox, compSamplesBox=argsList$compSamplesBox)
        }
        
        if(button$name == "geneSetsButton")
        {
            argsList$geneSetsTxt$setText(fileName)
        }
    }
    dialog$destroy()    
} 

############################################################################################################
# argsList= parentWindow, genesText, tabsTool, jobNameText
# gtLinker: 
# david: jobNameText, comboIdType, annotsArea, emailText, frameEmail, argsText
# gage:
# topgo:
# other:
# Common: tabsSteps, feaResultsText, comboFeaTool, expressionText
submitQuery <- function(button, argsList) # button: Not used, but required for signalConnect.
{     
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    activeTool <- rownames(FEA_tools)[which(FEA_tools[,"ID"]==RGtk2::gtkNotebookGetCurrentPage(argsList$tabFEA$commonFields$tabsTools))]
    geneList <- processGenesField(argsList$tabFEA$commonFields$genesText) 
    eset <- argsList$tabFEA$commonFields$esetTxt$getText()
    
    jobName <- argsList$tabFEA$commonFields$jobNameTex$getText() 
    if(jobName=="") jobName <- NULL    
    
    if(length(geneList)==0 && eset=="" )
    {
        # Dialog
        dialog <- RGtk2::gtkDialogNewWithButtons("Enter genes", argsList$parentWindow,
                                          c("modal", "destroy-with-parent"), 
                                          "gtk-ok", RGtk2::GtkResponseType["accept"],
                                          show=TRUE)
        if(activeTool=="gage") 
        {
            msgTxt <- "Please enter an expression set and select the samples to compare."
            
        } else 
        {
            msgTxt <- "Please enter a gene list."
        }
        dialog[["vbox"]]$add(RGtk2::gtkLabel(msgTxt))
        response <- dialog$run() 
        dialog$destroy()
        
    }else
    {
        msgText <- NULL
        statusTxt <- NULL
        if(activeTool == "Imported text file")
        {
            response <- RGtk2::GtkResponseType["accept"]
        }else
        {
            # Dialog    
            dialog <- RGtk2::gtkDialogNewWithButtons("Submit?", argsList$parentWindow,
                                             c("modal", "destroy-with-parent"), 
                                             "gtk-ok", RGtk2::GtkResponseType["accept"], 
                                             "gtk-cancel", RGtk2::GtkResponseType["reject"],
                                             show=FALSE)
            if(activeTool %in% c("GeneTerm Linker")) dialog[["vbox"]]$add(RGtk2::gtkLabel(paste("Submit query to ", activeTool, " server?", sep=""))) #"DAVID"
               if(activeTool %in% c("topGO", "gage")) dialog[["vbox"]]$add(RGtk2::gtkLabel("Perform analysis?"))
            response <- dialog$run() # showAll()
            dialog$destroy()        
        }    
        
        if (response == RGtk2::GtkResponseType["accept"]) 
        {    
            queryReply <- NULL
            argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Waiting for FEA results...")
            ##### GeneTerm Linker
            if(activeTool =="GeneTerm Linker")
            {
                # Check libraries (ask install through dialog)
                if(!loadInstPkg("RCurl", parentWindow=argsList$parentWindow)) stop("Package 'RCurl' is required to get the FEA results from GeneTerm Linker server.")
                
                queryArgs <- argsList$tabFEA$toolTabs$tabGTL$queryArgs    
                
                ## Get arguments queryArgs <- tabFEA$tabDavid$queryArgs 
                org <- names(queryArgs$GtL_org)[which(queryArgs$GtL_org==queryArgs$comboOrg$getActive())]
                anots <- queryArgs$GtL_annots[sapply(queryArgs$GtL_annots, function(anot) queryArgs$checkAnnotsGtL[[anot]]$active)]
                minSupport <- as.numeric(names(queryArgs$GtL_minSupport)[which(queryArgs$GtL_minSupport==queryArgs$comboMinSupport$getActive())]) 
                serverWeb <- queryArgs$serverWebText$getText()
                serverWS <- paste(serverWeb, ":8182", sep="")
            
                # Submit query
                queryReply <- -1
                queryReply <- tryCatch( 
                {
                    tmpReply <- fea_gtLinker(geneList=geneList, organism=org, annotations=anots, minSupport=minSupport, serverWS = serverWS)
                    tmpReply
                }, error = function(e) 
                {
                    msgText <<- paste("Error in query:\n " ,e, sep="")
                })
                
                if(is.null(msgText)) # No error
                {
                    if(queryReply==-1) msgText <- "There has been an error with the query. \nCheck the parameters and try again.\n"
                    if(queryReply!=-1)
                    {
                        msgText <- paste("Query submited.\nYour job ID is ",queryReply, ". \nIt will be ready in a few minutes.",sep="")
                        if(is.null(jobName)) jobName <- paste(queryReply, "_gtLinker",sep="")
                        statusTxt <- queryReply
                        argsList$serverWebReportText$setText(serverWeb)
                    }
                }
            }else
            {
                if(activeTool == "Imported text file")
                {
                     queryReply <- formatResultsFile(c(argsList$tabFEA$toolTabs$tabOther$fields, list(jobNameText=argsList$tabFEA$commonFields$jobNameText)))
                     msgText <- paste("The file has been formated to use with FGNet \n and saved as '", queryReply$fileName, "'\nYou can now build the functional network.",sep="")
                }
                # if(activeTool =="DAVID")
                # {                    
                #     queryArgs <- argsList$tabFEA$toolTabs$tabDavid$queryArgs    
                # 
                #     ## Get arguments
                #     geneIdType <- queryArgs$davidVars$DAV_GeneIds[queryArgs$comboIdType$getActive()+1]
                # 
                #     annotations <- unlist(sapply(queryArgs$annotsArea$getChildren(), function(anot) if(anot$active) return(anot$label)))
                #     email <- queryArgs$emailText$getText()
                #     if(email=="" || !RGtk2::gtkWidgetGetSensitive(queryArgs$frameEmail)) email <- NULL
                #     
                #     # Check libraries (ask install through dialog)
                #     if(!is.null(email))
                #     {
                #         if(!loadInstPkg("RDAVIDWebService", parentWindow=argsList$parentWindow)) 
                #             stop("Package RDAVIDWebService is required to query DAVID through the webserver. Install the package or set email=NULL to query DAVID through the web API.")
                #     }else{
                #         if(!loadInstPkg("RCurl", parentWindow=argsList$parentWindow)) stop("Package 'RCurl' is required to query DAVID throught the web API. Install it or provide an email to query DAVID through the Web Service.")
                #     }
                #         
                #     argsWS <- eval(parse(text=paste("c(",queryArgs$argsText$getText(), ")", sep="")))
                #     
                #     # Submit query
                #     queryReply <- tryCatch( 
                #     {                        
                #         tmpReply <- fea_david(geneList=geneList, geneIdType=geneIdType, annotations=annotations, email=email, argsWS=argsWS, jobName=jobName)
                #         tmpReply # (invisible returns)
                #     }, error = function(e) 
                #     {
                #         msgText <<- paste("Error in query:\n " ,e, sep="")
                #         NULL
                #     })
                # }
                    
                if(activeTool =="topGO")
                { 
                    # Check libraries (ask install through dialog)
                    if(!loadInstPkg("topGO", parentWindow=argsList$parentWindow)) stop("Package topGO is not available.")
                    if(!loadInstPkg("GO.db", parentWindow=argsList$parentWindow)) stop("Package GO.db is for fea_topGO.")
                    
                    queryArgs <- argsList$tabFEA$toolTabs$tabTopGo$queryArgs
                    
                    ## Get arguments
                    geneIdType <- queryArgs$geneIDsTopGO[queryArgs$comboGeneIdTypeTopGo$getActive()+1]
                    annotations <- names(queryArgs$topGO_annots)[sapply(queryArgs$topGO_annots, function(anot) queryArgs$checkAnnotsTopGo[[anot]]$active)]
                    organism <- names(queryArgs$allOrgs[["ID"]])[which(queryArgs$allOrgs[["ID"]]==queryArgs$comboOrgTopGo$getActive())]
                    evidence  <- rownames(queryArgs$GOEvidenceCodes)[sapply(rownames(queryArgs$GOEvidenceCodes), function(evid) queryArgs$checkEvidenceTopGo[[evid]]$active)]
                    if(length(evidence)==0) evidence <- NULL
                    
                    geneID2GO <- NULL           # geneID2GO <- buildDatabases(orgPackage, geneIdType)$geneID2GO
                    refPackage <- queryArgs$refPackageTxt$getText()
                    if(refPackage=="") refPackage <- NULL
                    genesUniverse <- NULL       # genesUniverse <- refList(refPackage, geneIdType)
                    nodeSize <- as.numeric(queryArgs$nodeSizeTxt$getText())
                    if(is.na(nodeSize)) nodeSize <- NULL
                    pValThr <- as.numeric(queryArgs$pValThrTxt$getText())
                    if(is.na(pValThr)) pValThr <- NULL
                    
                    # Submit query
                    queryReply <- tryCatch( 
                    {   
                        tmpReply <- fea_topGO(geneList=geneList, geneIdType=geneIdType, organism=organism, annotations=annotations, evidence=evidence, genesUniverse=genesUniverse, refPackage=refPackage, geneID2GO=geneID2GO, nodeSize=nodeSize, pValThr=pValThr, testStat=NULL)
                        tmpReply # (invisible returns)
                    }, error = function(e) 
                    {
                        msgText <<- paste("Error in query:\n " ,e, sep="")
                        NULL
                    })
                }    # End topGO
                
                
                if(activeTool =="gage")
                { 
                    # Check libraries (ask install through dialog)
                    if(!loadInstPkg("gage", parentWindow=argsList$parentWindow)) stop("Package 'gage' is not available.")
                    queryArgs <- argsList$tabFEA$toolTabs$tabGage$queryArgs    
                    
                    ## Get arguments
                    eset <- eval(as.name(load(eset)))
                    
                    refSamples <- unlist(sapply(argsList$tabFEA$commonFields$refSamplesBox$getChildren(), function(samp) if(samp$active) return(samp$label)))
                    compSamples <- unlist(sapply(argsList$tabFEA$commonFields$compSamplesBox$getChildren(), function(samp) if(samp$active) return(samp$label)))
                    
                    geneIdType <- queryArgs$geneIDsGage[queryArgs$comboGeneIdTypeGage$getActive()+1]
                    organism <- names(queryArgs$allOrgs[["ID"]])[which(queryArgs$allOrgs[["ID"]]==queryArgs$comboOrgGage$getActive())]
                    
                    annotations <- names(queryArgs$gage_annots)[sapply(queryArgs$gage_annots, function(anot) queryArgs$checkAnnotsGage[[anot]]$active)]
                    if(length(annotations)==0) annotations<-NULL
                    geneSets <- queryArgs$geneSetsTxt$getText()
                    if(geneSets=="") geneSets <- NULL
                    if(!is.null(geneSets))
                    {
                        # http://www.broadinstitute.org/gsea/msigdb/collections.jsp
                        if(grepl(".gmt", geneSets, fixed=TRUE))
                        {
                            geneSets <- readList(geneSets)
                        }else{
                            geneSets <- eval(as.name(load(geneSets)))    
                        }
                    }
                    
                    sameDirection <- queryArgs$sameDirectionCheck$active
                    onlyEssentialTerms <- queryArgs$onlyEssentialTermsCheck$active
                    compareType <- queryArgs$compareTypeTxt$getText()
                    
                    # Submit query
                    queryReply <- tryCatch( 
                    {                                            
                        tmpReply <- fea_gage(eset=eset, refSamples=refSamples, compSamples=compSamples, geneIdType=geneIdType, organism=organism, annotations=annotations, geneSets=geneSets, sameDirection=sameDirection, onlyEssentialTerms=onlyEssentialTerms, compareType=compareType, jobName=jobName)
                        tmpReply # (invisible returns)
                            
                    }, error = function(e) 
                    {
                        msgText <<- paste("Error in query:\n " ,e, sep="")
                        NULL
                    })
                }    # End gage
                
                if(is.null(msgText)) # No error
                {
                    if(is.null(jobName)) jobName <- getJobName(queryReply$fileName)
                                        
                    if(!is.null(queryReply)) 
                    {
                        msgText <- paste("The results from the Functional Analysis have been saved in folder '",jobName, "'.\nYou can now build the functional network.",sep="")
                       
                    }
                    if(is.null(queryReply)) 
                    {
                        # In case there isn't an specific msg from the tool...
                        msgText <- "There has been an error with the query. \nCheck the parameters and try again."
                        queryReply <- NULL
                    }
                }
                statusTxt <- jobName
            }# End queries
            argsList$statusbar$push(argsList$statusbar$getContextId("info"), paste("FEA ready: ", statusTxt, sep=""))
        }
    
        if(!is.null(msgText))
        {
            ############
            # Comun
            # queryReply dialog:
            dialog <- RGtk2::gtkDialogNewWithButtons("Query submited", argsList$parentWindow,
                                              c("modal", "destroy-with-parent"), 
                                              "gtk-ok", RGtk2::GtkResponseType["accept"], 
                                              show=FALSE)
            dialog[["vbox"]]$add(RGtk2::gtkLabel(msgText))
            response <- dialog$run() # showAll()
            dialog$destroy()                

            if(is.numeric(queryReply) && queryReply==-1) queryReply <- NULL # Error 
            
            if(!is.null(queryReply)) # 
            {                
                if(is.list(queryReply)) # (Not GTLinker)
                {
                    if(is.null(jobName)) jobName <- getJobName(queryReply$fileName)
                    feaResults <- queryReply
                    queryReply <- paste(jobName, .Platform$file.sep, jobName, "_feaResults.RData", sep="") 
                    save(feaResults, file=queryReply)
                }
                

                # Go to Network tab and pass values
                RGtk2::gtkNotebookSetCurrentPage(argsList$tabsSteps, argsList$tabID["Network"])
                argsList$feaResultsText$setText(queryReply)
                argsList$comboFeaTool$setActive(FEA_tools[activeTool,"ID"])
                
                if("genesFC" %in% names(feaResults)) 
                {
                    genesFC <- sort(names(feaResults$genesFC))                    

                    expressionBuffer <- RGtk2::gtkTextBufferNew()
                    RGtk2::gtkTextBufferSetText(expressionBuffer,  paste(genesFC,round(feaResults$genesFC[genesFC], digits=3), sep="\t", collapse="\n"))
                    RGtk2::gtkTextViewSetBuffer(argsList$expressionText, expressionBuffer)                    
                }
            }
            msgText <- NULL
        }
    }
}   


############################################################################################################
## Get FEA results
readFEAresults <- function(feaResultsText, comboFeaTool, expressionText, serverWebReportText, alreadyDownloadedCheck, parentWindow)
{
    gtLinkerWsPort <- ":8182"
    
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    fileName <- feaResultsText$getText()
    
    if(grepl(".RData",fileName, fixed=TRUE))
    {
        feaResults <- eval(as.name(load(fileName)))
        if("genesFC" %in% names(feaResults)) 
        {
            genesFC <- sort(names(feaResults$genesFC))                    
            
            expressionBuffer <- RGtk2::gtkTextBufferNew()
            RGtk2::gtkTextBufferSetText(expressionBuffer,  paste(genesFC,round(feaResults$genesFC[genesFC], digits=3), sep="\t", collapse="\n"))
            RGtk2::gtkTextViewSetBuffer(expressionText, expressionBuffer)                    
        }
    }else{
        reportTool <- comboFeaTool$getActive()
        reportTool <- rownames(FEA_tools)[which(as.numeric(FEA_tools[,"ID"])==reportTool)]
        
        if(reportTool == "GeneTerm Linker")
        {
            jobID <- fileName
            alreadyDownloaded <- alreadyDownloadedCheck$active
            serverWeb <- serverWebReportText$getText()
            serverWS <- paste(serverWebReportText$getText(), gtLinkerWsPort, sep="")
            # setwd("..") # GtLinker moves into the folder
            msgText <- NULL
            feaResults <- tryCatch({
                tmpResults <- suppressWarnings(fea_gtLinker_getResults(jobID=jobID, jobName=NULL, alreadyDownloaded=alreadyDownloaded,serverWeb=serverWeb,serverWS=serverWS, keepTrying=FALSE))
                tmpResults
            }, error = function(e) 
            {
                msgText <<- paste("Error in query:\n " ,e, sep="")
                NULL
            })
            if(is.null(feaResults))
            {
                msgText <- "The analysis has not finished yet or the jobID does not exist."
            }
            
            if(!is.null(msgText))
            {
                dialog <- RGtk2::gtkDialogNewWithButtons("FEA not ready", parentWindow,
                                                  c("modal", "destroy-with-parent"), 
                                                  "gtk-ok", RGtk2::GtkResponseType["accept"], 
                                                  show=FALSE)
                dialog[["vbox"]]$add(RGtk2::gtkLabel(msgText))
                response <- dialog$run() 
                dialog$destroy() 
                #stop()
            }else{
                # Save and replace jobID by the downloaded file:
                jobName <- getJobName(feaResults$fileName)
                fileName <- paste(jobName,"_feaResults.RData", sep="")
                fileName <- paste(getwd(), jobName, fileName, sep=.Platform$file.sep)
                save(feaResults, file=fileName)
                feaResultsText$setText(fileName) 
                alreadyDownloadedCheck$active <- TRUE
                
                
                # Add mgs to combo in subnetwork tab
                # subnw:
                #                 clSelectClusterFrame$label <- "Select metagroup"
                #                 clustersId <<- setNames(as.numeric(rownames(feaResults$metagroups))-1, rownames(feaResults$metagroups))
            }
            
        }else
        {
            feaResults <- readGeneTermSets(fileName, tool=reportTool,simplifyGage=FALSE) # geneLabels=NULL?
            
            # subnw:
            #clSelectClusterFrame$label <- paste("Select ", tolower(FEA_tools[reportTool,"GroupType"]), sep="")
            #             if(reportTool %in% c("DAVID","gage", "Imported text file"))
            #             {
            #                 # Add cl to combo in subnetwork tab
            #                 clustersId <<- setNames(0:(length(feaResults$clusters$Cluster)-1), feaResults$clusters$Cluster)
            #             } 
            #             if(reportTool %in% c("topGO"))
            #             {
            #                 # Add cl to combo in subnetwork tab
            #                 clustersId <<- setNames(1:nrow(feaResults)-1, 1:nrow(feaResults))
            #             }
        }
    }
    return(feaResults)
}

