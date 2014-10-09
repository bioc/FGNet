################################################################################
### Aux functions
# selectNetworkType
# generateNetwork

# argsList= comboFunNwIntNw, funNwIntNwID, plotTitleText, weightedCheck, keepAllNodesCheck, comboPlotOutput,plotOutputID
selectNetworkType <- function(aux, argsList)
{
  selNwTypeID<- which(argsList$funNwIntNwID==argsList$comboFunNwIntNw$getActive())
  selNwTypeTxt <- names(argsList$funNwIntNwText)[selNwTypeID]
  
  if(names(selNwTypeID) %in% c("default","bipartite"))
  {
    argsList$plotTitleText$setText(selNwTypeTxt)
    gtkWidgetSetSensitive(argsList$plotTitleText, TRUE)
    gtkWidgetSetSensitive(argsList$comboPlotOutput, TRUE)
    
    if(names(selNwTypeID) == "default")
    {
      gtkWidgetSetSensitive(argsList$weightedCheck, TRUE)
      gtkWidgetSetSensitive(argsList$keepAllNodesCheck, FALSE)
    }
    if(names(selNwTypeID) == "bipartite")
    {
      gtkWidgetSetSensitive(argsList$weightedCheck, FALSE)
      gtkWidgetSetSensitive(argsList$keepAllNodesCheck, TRUE)
    }
    
  }
  if(names(selNwTypeID) == "both")
  {
    argsList$comboPlotOutput$setActive(argsList$plotOutputID["dynamic"])
    gtkWidgetSetSensitive(argsList$comboPlotOutput, FALSE)
    argsList$plotTitleText$setText("")
    gtkWidgetSetSensitive(argsList$plotTitleText, FALSE)
    
    gtkWidgetSetSensitive(argsList$weightedCheck, TRUE)
    gtkWidgetSetSensitive(argsList$keepAllNodesCheck, TRUE)
  }
}


# argsList= 
# parentWindow, statusbar
# commonFields: feaResultsText, comboFeaTool,  thresholdAttributeNetworkText, thresholdValueNetworkText, thresholdOperatorNetworkText, expressionText, comboExpressionPlot
# all tabPlotNetwork fields (onlyNwFields)

generateNetwork <- function(button, argsList) # button: Not used, but required for signalConnect.
{ 
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    ###################################################################
    # Start processing
    inputFile <- argsList$commonFields$feaResultsText$getText()
    reportTool <- argsList$commonFields$comboFeaTool$getActive()
    reportTool <- rownames(FEA_tools)[which(as.numeric(FEA_tools[,"ID"])==reportTool)]
    
    # Recargar siempre? donde poner para q solo cuando cambia?
    feaResults <- readFEAresults(feaResultsText=argsList$commonFields$feaResultsText,comboFeaTool=argsList$commonFields$comboFeaTool, expressionText=argsList$commonFields$expressionText,
                                    serverWebReportText=argsList$commonFields$serverWebReportText,alreadyDownloadedCheck=argsList$commonFields$alreadyDownloadedCheck, parentWindow=argsList$parentWindow)
    
    # inputFile is GTLinker jobID?
    if(!is.na(suppressWarnings(as.numeric(inputFile)))) inputFile <- feaResults$fileName
        
    if(is.null(feaResults) || (reportTool==-1 || inputFile==""))
    {
        # Dialog
        dialog <- gtkDialogNewWithButtons("FEA results required", argsList$parentWindow,
                                          c("modal", "destroy-with-parent"), 
                                          "gtk-ok", GtkResponseType["accept"],
                                          show=TRUE)
        dialog[["vbox"]]$add(gtkLabel("Please provide a file with the FEA results to continue."))
        response <- dialog$run() 
        dialog$destroy()
        
    }else
    {        
        ###################################################################
        # Get arguments
        geneExpr <- processExpressionField(argsList$commonFields$expressionText)
        plotExpression <- names(argsList$commonFields$exprPlotOptionsID)[which(argsList$commonFields$exprPlotOptionsID==argsList$commonFields$comboExpressionPlot$getActive())]
        ###################################################################
        # Go to the given file directory        
        folder <- strsplit(inputFile, .Platform$file.sep, fixed=TRUE)[[1]]
        folder <- sub(folder[length(folder)],"",inputFile)
                      
        # Create folder
        currWD <- getwd()
        if(folder!="")
        {
            if((!file.exists(file.path(folder))))
            {
                dir.create(file.path(folder))        
            }
            setwd(folder)
        }
        
        ####################################
        #### Plot networks
        retFNw <- NULL
        if(!is.null(feaResults))
        {
            
          argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Generating network...")
          # Subnetwork
          if(button$label %in% c("nw1","nw2","nw3","nw5","nw6","nw7","nw8"))
          {
#               # Get selected cluster
#               selCluster=""
#               if(txtSelectCluster$visible) selCluster <- txtSelectCluster$getText()
#               if(comboSelectCluster$visible) selCluster <- names(clustersId)[which(clustersId==comboSelectCluster$getActive())]
#               if(selCluster=="") selCluster <- clustersId[1]
#               txtSelectCluster$setText("")
#               
#               # Fill combo          
#               comboSelectCluster$getModel()$clear()
#               for (i in names(clustersId)) gtkComboBoxAppendText(comboSelectCluster, i)
#               comboSelectCluster$setActive(clustersId[selCluster])
#               gtkWidgetSetSensitive(comboSelectCluster, TRUE)
#               comboSelectCluster$Show()
#               gtkWidgetSetSensitive(txtSelectCluster, FALSE)
#               txtSelectCluster$Hide()  
#               
#             subNwType <- button$label
#             plotOutput <- names(argsList$onlyNwFields$plotOutputID)[which(argsList$onlyNwFields$plotOutputID==argsList$onlyNwFields$comboClplotOutput$getActive())] # "dinamic","static"
# 
              ## add trycatch
#             retFNw <- plotSubNetwork(feaResults, selCluster=selCluster, plotOutput=subNwType, plotOutput=plotOutput, returnGraph=saveCheck$active, geneExpr=geneExpr,plotExpression=plotExpression)
#             
          }else # FunctionalNetwork 
          {
            jobName <- getJobName(feaResults$fileName) # for save
            ##############################
            # Get arguments
            # toMatrix():
            key <- names(argsList$onlyNwFields$keysID)[which(argsList$onlyNwFields$keysID==argsList$onlyNwFields$comboKey$getActive())]
            sepChar <- NULL
            clusterColumn <- NULL
            filterAttribute <- argsList$commonFields$thresholdAttributeNetworkText$getText()  # Only david and gTLinker
                if(filterAttribute=="") filterAttribute <- NULL
            filterThreshold <- argsList$commonFields$thresholdValueNetworkText$getText()
                if(filterThreshold=="") filterThreshold <- 0
            filterOperator <- argsList$commonFields$thresholdOperatorNetworkText$getText() 
                if(filterOperator=="") filterOperator <- NULL
            removeFiltered <- NULL
            
            # functionalNetwork():
            weighted <- argsList$onlyNwFields$weightedCheck$active
            keepAllNodes <- argsList$onlyNwFields$keepAllNodesCheck$active
            plotAllMg <- keepAllNodes
            plotTitle <- argsList$onlyNwFields$plotTitleText$getText()
            legendPrefix <- argsList$onlyNwFields$legendPrefixText$getText()
            legendText <- NULL
            if(!argsList$onlyNwFields$showLegendCheck$active) legendText <- FALSE
            vSize <- as.numeric(argsList$onlyNwFields$vSizeText$getText())
            vLabelCex <- as.numeric(argsList$onlyNwFields$vLabelCexText$getText())
            bgTransparency <- as.numeric(argsList$onlyNwFields$bgTransparencyText$getText())
            eColor <- argsList$onlyNwFields$eColorText$getText()
            nwType <- names(argsList$onlyNwFields$funNwIntNwID)[which(argsList$onlyNwFields$funNwIntNwID==argsList$onlyNwFields$comboFunNwIntNw$getActive())]  # "functionalNetwork","intersectionNetwork","both"
            if(nwType=="both") nwType <- names(argsList$onlyNwFields$funNwIntNwID)[1:2]
            
            plotOutput <- names(argsList$onlyNwFields$plotOutputID)[which(argsList$onlyNwFields$plotOutputID==argsList$onlyNwFields$comboPlotOutput$getActive())] # "dinamic","static"
            
            #### PLOT
            # Matrices
            matr <- fea2incidMat(feaResults=feaResults, key=key, sepChar=sepChar, clusterColumn=clusterColumn, filterAttribute=filterAttribute, filterThreshold=filterThreshold,filterOperator=filterOperator, removeFilteredGtl=removeFiltered)
            
            # Plot
            retFNw <- tryCatch( 
                            {                        
                                tmpReply <- functionalNetwork(matr, plotOutput=plotOutput, plotType=nwType, 
                                                              plotTitle=plotTitle, legendPrefix=legendPrefix, legendText=legendText, 
                                                              plotExpression=plotExpression,geneExpr=geneExpr, 
                                                              vSize=vSize, vLabelCex=vLabelCex, eColor=eColor, 
                                                              weighted=weighted,  bgTransparency=bgTransparency, 
                                                              plotTitleSub=NULL, keepAllNodes=keepAllNodes, plotAllMg=plotAllMg)  # MISSING!!!  
                                argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Network ready.")
                                tmpReply # (invisible returns)
                            }, error = function(e) 
                            {
                                argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Error in network. See R console for details.")
                                stop(e)
                            })

          }
          
          ####################################
          # Save return
          if(argsList$onlyNwFields$saveCheck$active)
          {
              if(jobName=="") jobName <- "network"
              jobName <- paste(jobName,"_network.RData",sep="")
              save(retFNw, file=jobName)
              message(paste("Network saved as ",jobName, sep=""))
          }
        }
        setwd(currWD)
    }    
} 

#################################################################################
### Fills the Only network tab
#################################################################################

tabPlotNetwork_fill <- function(mainWindow, statusbar, commonFields)
{
    tabPlotNetwork <- gtkVBox(FALSE,3)
    ### plotOptionsFrame
    # Not added: sepChar, vLayout=NULL, keepColors=TRUE
    plotOptionsExpander <- gtkExpanderNew("More options")
    fgnPlotOptionsBox <- gtkVBox(FALSE,3)
    
    fgnPlotOptionsL1Box <- gtkHBox(FALSE,3)
    fgnPlotOptionsL2Box <- gtkHBox(FALSE,3)
    fgnPlotOptionsL1bBox <- gtkHBox(FALSE,3)
    fgnPlotOptionsL3Box <- gtkHBox(FALSE,3)
    fgnPlotOptionsL4Box <- gtkHBox(FALSE,3)
    fgnPlotOptionsL5Box <- gtkHBox(FALSE,3)
    
    # key="Genes"
    keysID <- setNames(0:1, c("Genes","Terms"))
    comboKeyFrame <- gtkFrame("Network type")
    gtkFrameSetShadowType(comboKeyFrame, GtkShadowType["none"])
    comboKey <- gtkComboBoxNewText()
    for (i in names(keysID)) gtkComboBoxAppendText(comboKey, i)
    comboKey$setActive(keysID["Genes"])
    comboKeyFrame$add(comboKey)
    fgnPlotOptionsL1Box$packStart(comboKeyFrame, expand=FALSE)
    
    # plotFunctionalNetwork
    funNwIntNwID <- setNames(0:2, c("default","bipartite","both"))
    funNwIntNwText <- setNames(0:2, c("Functional Network (default)","Intersection Network (bipartite)","Both"))
    funNwIntNwFrame <- gtkFrame("")
    gtkFrameSetShadowType(funNwIntNwFrame, GtkShadowType["none"])
    comboFunNwIntNw <- gtkComboBoxNewText()
    for (i in names(funNwIntNwText)) gtkComboBoxAppendText(comboFunNwIntNw, i)
    comboFunNwIntNw$setActive(funNwIntNwID["default"])
    funNwIntNwFrame$add(comboFunNwIntNw)    
    fgnPlotOptionsL1Box$packStart(funNwIntNwFrame, expand=FALSE)
    
    # plotOutput="static"
    plotOutputID <- setNames(0:1, c("dynamic","static"))
    plotOutputFrame <- gtkFrame("Output")
    gtkFrameSetShadowType(plotOutputFrame, GtkShadowType["none"])
    comboPlotOutput <- gtkComboBoxNewText()
    for (i in capitalize(names(plotOutputID)))gtkComboBoxAppendText(comboPlotOutput, i)
    comboPlotOutput$setActive(plotOutputID["static"])
    plotOutputFrame$add(comboPlotOutput)
    fgnPlotOptionsL2Box$packStart(plotOutputFrame, expand=FALSE)
    
    # Checkbox
    fgnPlotOptionsHBox <- gtkHBox(FALSE,3)
    # save
    saveCheck<- gtkCheckButton("Save (.RData)")
    saveCheck$"tooltip-text" <- "Save igraph & matrices"
    fgnPlotOptionsL2Box$packStart(saveCheck, expand=FALSE)
    # weighted=FALSE
    weightedCheck<- gtkCheckButton("Weighted")
    fgnPlotOptionsHBox$packStart(weightedCheck, expand=FALSE)
    # keepAllNodes=FALSE
    keepAllNodesCheck<- gtkCheckButton("keepAllNodes")
    gtkWidgetSetSensitive(keepAllNodesCheck, FALSE)
    
    
    fgnPlotOptionsHBox$packStart(keepAllNodesCheck, expand=FALSE)
    fgnPlotOptionsL1bBox$packStart(fgnPlotOptionsHBox, expand=FALSE)
    
    
    
    # plotTitle="Functional Network"
    plotTitleFrame <- gtkFrame("Title")
    gtkFrameSetShadowType(plotTitleFrame, GtkShadowType["none"])
    plotTitleText <- gtkEntryNew()
    plotTitleText$setWidthChars(15)
    plotTitleText$setText("Functional Network")
    plotTitleFrame$add(plotTitleText)
    fgnPlotOptionsL3Box$packStart(plotTitleFrame, expand=TRUE)
    
    # legendPrefix=NULL
    legendPrefixFrame <- gtkFrame("Legend prefix")
    legendPrefixBox <- gtkHBox(FALSE,3)
    gtkFrameSetShadowType(legendPrefixFrame, GtkShadowType["none"])
    legendPrefixText <- gtkEntryNew()
    legendPrefixText$setText("Cl ")
    legendPrefixText$setWidthChars(7)
    legendPrefixBox$packStart(legendPrefixText, expand=FALSE)
    
    showLegendCheck<- gtkCheckButton("Show legend")
    showLegendCheck$active <- TRUE
    legendPrefixBox$packStart(showLegendCheck, expand=FALSE)
    
    legendPrefixFrame$add(legendPrefixBox)
    fgnPlotOptionsL3Box$packStart(legendPrefixFrame, expand=FALSE)
    
    # legendText=NULL
    
    # vSize=12
    vSizeFrame <- gtkFrame("vSize")
    gtkFrameSetShadowType(vSizeFrame, GtkShadowType["none"])
    vSizeText <- gtkEntryNew()
    vSizeText$setWidthChars(5)
    vSizeText$setText("12")
    vSizeFrame$add(vSizeText)
    fgnPlotOptionsL4Box$packStart(vSizeFrame, expand=TRUE)
    
    # vLabelCex=0.75 
    vLabelCexFrame <- gtkFrame("LabelCex")
    gtkFrameSetShadowType(vLabelCexFrame, GtkShadowType["none"])
    vLabelCexText <- gtkEntryNew()
    vLabelCexText$setWidthChars(5)
    vLabelCexText$setText("0.75")
    vLabelCexFrame$add(vLabelCexText)
    fgnPlotOptionsL4Box$packStart(vLabelCexFrame, expand=TRUE)
    
    # bgTransparency=0.4
    bgTransparencyFrame <- gtkFrame("bg")
    gtkFrameSetShadowType(bgTransparencyFrame, GtkShadowType["none"])
    bgTransparencyText <- gtkEntryNew()
    bgTransparencyText$"tooltip-text" <- "Background transparency"
    bgTransparencyText$setWidthChars(5)
    bgTransparencyText$setText("0.4")
    bgTransparencyFrame$add(bgTransparencyText)
    fgnPlotOptionsL4Box$packStart(bgTransparencyFrame, expand=TRUE)
    
    # eColor="#323232"
    eColorFrame <- gtkFrame("eColor")
    gtkFrameSetShadowType(eColorFrame, GtkShadowType["none"])
    eColorText <- gtkEntryNew()
    bgTransparencyText$setWidthChars(5)
    eColorText$setText("#323232")
    eColorFrame$add(eColorText)
    fgnPlotOptionsL4Box$packStart(eColorFrame, expand=FALSE)
    
    # Tab end
    fgnPlotOptionsBox$packStart(fgnPlotOptionsL3Box, expand=FALSE,10)
    fgnPlotOptionsBox$packStart(fgnPlotOptionsL4Box, expand=FALSE,10)
    fgnPlotOptionsBox$packStart(fgnPlotOptionsL5Box, expand=FALSE,10)
    plotOptionsExpander$add(fgnPlotOptionsBox)
    
    tabPlotNetwork$packStart(fgnPlotOptionsL1Box, expand=FALSE,10)
    tabPlotNetwork$packStart(fgnPlotOptionsL1bBox, expand=FALSE,10)
    tabPlotNetwork$packStart(fgnPlotOptionsL2Box, expand=FALSE,10)
    tabPlotNetwork$packStart(plotOptionsExpander, expand=FALSE,10)
    
    buttonNetwork <- gtkButton("Create")
    tabPlotNetwork$packStart(buttonNetwork, expand = FALSE, FALSE, 10) 
    
    # Make list with fields:
    fields <- list(comboKey=comboKey,keysID=keysID, comboFunNwIntNw=comboFunNwIntNw, funNwIntNwID=funNwIntNwID, funNwIntNwText=funNwIntNwText, weightedCheck=weightedCheck, keepAllNodesCheck=keepAllNodesCheck, comboPlotOutput=comboPlotOutput,plotOutputID=plotOutputID, saveCheck=saveCheck,
    plotTitleText=plotTitleText, legendPrefixText=legendPrefixText, showLegendCheck=showLegendCheck, 
    vSizeText=vSizeText, vLabelCexText=vLabelCexText, bgTransparencyText=bgTransparencyText, eColorText=eColorText)
    
    
    # SignalConnects
    gSignalConnect(buttonNetwork, "clicked", generateNetwork, data=list(parentWindow=mainWindow, statusbar=statusbar, 
                                                            commonFields=commonFields, onlyNwFields=fields))
    
    gSignalConnect(comboFunNwIntNw, "changed", selectNetworkType, data=list(comboFunNwIntNw=comboFunNwIntNw,funNwIntNwID=funNwIntNwID, funNwIntNwText=funNwIntNwText, plotTitleText=plotTitleText, comboPlotOutput=comboPlotOutput,plotOutputID=plotOutputID,
                                                                            weightedCheck=weightedCheck, keepAllNodesCheck=keepAllNodesCheck))
    
    return(list(tabPlotNetwork=tabPlotNetwork, fields=fields))
}