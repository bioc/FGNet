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
    RGtk2::gtkWidgetSetSensitive(argsList$plotTitleText, TRUE)
    RGtk2::gtkWidgetSetSensitive(argsList$comboPlotOutput, TRUE)
    
    if(names(selNwTypeID) == "default")
    {
      RGtk2::gtkWidgetSetSensitive(argsList$weightedCheck, TRUE)
      RGtk2::gtkWidgetSetSensitive(argsList$keepAllNodesCheck, FALSE)
    }
    if(names(selNwTypeID) == "bipartite")
    {
      RGtk2::gtkWidgetSetSensitive(argsList$weightedCheck, FALSE)
      RGtk2::gtkWidgetSetSensitive(argsList$keepAllNodesCheck, TRUE)
    }
    
  }
  if(names(selNwTypeID) == "both")
  {
    argsList$comboPlotOutput$setActive(argsList$plotOutputID["dynamic"])
    RGtk2::gtkWidgetSetSensitive(argsList$comboPlotOutput, FALSE)
    argsList$plotTitleText$setText("")
    RGtk2::gtkWidgetSetSensitive(argsList$plotTitleText, FALSE)
    
    RGtk2::gtkWidgetSetSensitive(argsList$weightedCheck, TRUE)
    RGtk2::gtkWidgetSetSensitive(argsList$keepAllNodesCheck, TRUE)
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
        dialog <- RGtk2::gtkDialogNewWithButtons("FEA results required", argsList$parentWindow,
                                          c("modal", "destroy-with-parent"), 
                                          "gtk-ok", RGtk2::GtkResponseType["accept"],
                                          show=TRUE)
        dialog[["vbox"]]$add(RGtk2::gtkLabel("Please provide a file with the FEA results to continue."))
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
#               for (i in names(clustersId)) RGtk2::gtkComboBoxAppendText(comboSelectCluster, i)
#               comboSelectCluster$setActive(clustersId[selCluster])
#               RGtk2::gtkWidgetSetSensitive(comboSelectCluster, TRUE)
#               comboSelectCluster$Show()
#               RGtk2::gtkWidgetSetSensitive(txtSelectCluster, FALSE)
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
    tabPlotNetwork <- RGtk2::gtkVBox(FALSE,3)
    ### plotOptionsFrame
    # Not added: sepChar, vLayout=NULL, keepColors=TRUE
    plotOptionsExpander <- RGtk2::gtkExpanderNew("More options")
    fgnPlotOptionsBox <- RGtk2::gtkVBox(FALSE,3)
    
    fgnPlotOptionsL1Box <- RGtk2::gtkHBox(FALSE,3)
    fgnPlotOptionsL2Box <- RGtk2::gtkHBox(FALSE,3)
    fgnPlotOptionsL1bBox <- RGtk2::gtkHBox(FALSE,3)
    fgnPlotOptionsL3Box <- RGtk2::gtkHBox(FALSE,3)
    fgnPlotOptionsL4Box <- RGtk2::gtkHBox(FALSE,3)
    fgnPlotOptionsL5Box <- RGtk2::gtkHBox(FALSE,3)
    
    # key="Genes"
    keysID <- setNames(0:1, c("Genes","Terms"))
    comboKeyFrame <- RGtk2::gtkFrame("Network type")
    RGtk2::gtkFrameSetShadowType(comboKeyFrame, GtkShadowType["none"])
    comboKey <- RGtk2::gtkComboBoxNewText()
    for (i in names(keysID)) RGtk2::gtkComboBoxAppendText(comboKey, i)
    comboKey$setActive(keysID["Genes"])
    comboKey$"tooltip-text" <- "Network nodes"
    comboKeyFrame$add(comboKey)
    fgnPlotOptionsL1Box$packStart(comboKeyFrame, expand=FALSE)
    
    # plotFunctionalNetwork
    funNwIntNwID <- setNames(0:2, c("default","bipartite","both"))
    funNwIntNwText <- setNames(0:2, c("Functional Network (default)","Intersection Network (bipartite)","Both"))
    funNwIntNwFrame <- RGtk2::gtkFrame("")
    RGtk2::gtkFrameSetShadowType(funNwIntNwFrame, GtkShadowType["none"])
    comboFunNwIntNw <- RGtk2::gtkComboBoxNewText()
    for (i in names(funNwIntNwText)) RGtk2::gtkComboBoxAppendText(comboFunNwIntNw, i)
    comboFunNwIntNw$setActive(funNwIntNwID["default"])
    funNwIntNwFrame$add(comboFunNwIntNw)    
    fgnPlotOptionsL1Box$packStart(funNwIntNwFrame, expand=FALSE)
    
    # plotOutput="static"
    plotOutputID <- setNames(0:1, c("dynamic","static"))
    plotOutputFrame <- RGtk2::gtkFrame("Output")
    RGtk2::gtkFrameSetShadowType(plotOutputFrame, GtkShadowType["none"])
    comboPlotOutput <- RGtk2::gtkComboBoxNewText()
    for (i in capitalize(names(plotOutputID))) RGtk2::gtkComboBoxAppendText(comboPlotOutput, i)
    comboPlotOutput$setActive(plotOutputID["static"])
    comboPlotOutput$"tooltip-text" <- "static: standard R plot (R console) \ndynamic:interactive tkplot (metagroups background cannot be drawn)"
    
    plotOutputFrame$add(comboPlotOutput)
    fgnPlotOptionsL2Box$packStart(plotOutputFrame, expand=FALSE)
    
    # Checkbox
    fgnPlotOptionsHBox <- RGtk2::gtkHBox(FALSE,3)
    # save
    saveCheck<- RGtk2::gtkCheckButton("Save (.RData)")
    saveCheck$"tooltip-text" <- "Save igraph & matrices"
    fgnPlotOptionsL2Box$packStart(saveCheck, expand=FALSE)
    # weighted=FALSE
    weightedCheck<- RGtk2::gtkCheckButton("Weighted")
    weightedCheck$"tooltip-text" <- "Edge width based on the number of shared gene-term sets"    
    fgnPlotOptionsHBox$packStart(weightedCheck, expand=FALSE)
    # keepAllNodes=FALSE
    keepAllNodesCheck<- RGtk2::gtkCheckButton("keepAllNodes")
    keepAllNodesCheck$"tooltip-text" <- "Show all nodes, including those in only one cluster"
    RGtk2::gtkWidgetSetSensitive(keepAllNodesCheck, FALSE)
    
    
    fgnPlotOptionsHBox$packStart(keepAllNodesCheck, expand=FALSE)
    fgnPlotOptionsL1bBox$packStart(fgnPlotOptionsHBox, expand=FALSE)
    
    
    
    # plotTitle="Functional Network"
    plotTitleFrame <- RGtk2::gtkFrame("Title")
    RGtk2::gtkFrameSetShadowType(plotTitleFrame, GtkShadowType["none"])
    plotTitleText <- RGtk2::gtkEntryNew()
    plotTitleText$setWidthChars(15)
    plotTitleText$setText("Functional Network")
    plotTitleFrame$add(plotTitleText)
    fgnPlotOptionsL3Box$packStart(plotTitleFrame, expand=TRUE)
    
    # legendPrefix=NULL
    legendPrefixFrame <- RGtk2::gtkFrame("Legend prefix")
    legendPrefixBox <- RGtk2::gtkHBox(FALSE,3)
    RGtk2::gtkFrameSetShadowType(legendPrefixFrame, GtkShadowType["none"])
    legendPrefixText <- RGtk2::gtkEntryNew()
    legendPrefixText$setText("Cl ")
    legendPrefixText$setWidthChars(7)
    legendPrefixBox$packStart(legendPrefixText, expand=FALSE)
    
    showLegendCheck<- RGtk2::gtkCheckButton("Show legend")
    showLegendCheck$active <- TRUE
    legendPrefixBox$packStart(showLegendCheck, expand=FALSE)
    
    legendPrefixFrame$add(legendPrefixBox)
    fgnPlotOptionsL3Box$packStart(legendPrefixFrame, expand=FALSE)
    
    # legendText=NULL
    
    # vSize=12
    vSizeFrame <- RGtk2::gtkFrame("vSize")
    RGtk2::gtkFrameSetShadowType(vSizeFrame, GtkShadowType["none"])
    vSizeText <- RGtk2::gtkEntryNew()
    vSizeText$setWidthChars(5)
    vSizeText$setText("12")
    vSizeText$"tooltip-text" <- "Vertex (node) size"
    vSizeFrame$add(vSizeText)
    fgnPlotOptionsL4Box$packStart(vSizeFrame, expand=TRUE)
    
    # vLabelCex=0.75 
    vLabelCexFrame <- RGtk2::gtkFrame("LabelCex")
    RGtk2::gtkFrameSetShadowType(vLabelCexFrame, GtkShadowType["none"])
    vLabelCexText <- RGtk2::gtkEntryNew()
    vLabelCexText$setWidthChars(5)
    vLabelCexText$setText("0.75")
    vLabelCexText$"tooltip-text" <- "Node label size"
    vLabelCexFrame$add(vLabelCexText)
    fgnPlotOptionsL4Box$packStart(vLabelCexFrame, expand=TRUE)
    
    # bgTransparency=0.4
    bgTransparencyFrame <- RGtk2::gtkFrame("bgTransp.")
    RGtk2::gtkFrameSetShadowType(bgTransparencyFrame, GtkShadowType["none"])
    bgTransparencyText <- RGtk2::gtkEntryNew()
    bgTransparencyText$"tooltip-text" <- "Background transparency"
    bgTransparencyText$setWidthChars(5)
    bgTransparencyText$setText("0.4")
    bgTransparencyFrame$add(bgTransparencyText)
    fgnPlotOptionsL4Box$packStart(bgTransparencyFrame, expand=TRUE)
    
    # eColor="#323232"
    eColorFrame <- RGtk2::gtkFrame("eColor")
    RGtk2::gtkFrameSetShadowType(eColorFrame, GtkShadowType["none"])
    eColorText <- RGtk2::gtkEntryNew()
    bgTransparencyText$setWidthChars(5)
    eColorText$setText("#323232")
    eColorText$"tooltip-text" <- "Edge color"
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
    
    buttonNetwork <- RGtk2::gtkButton("Create")
    tabPlotNetwork$packStart(buttonNetwork, expand = FALSE, FALSE, 10) 
    
    # Make list with fields:
    fields <- list(comboKey=comboKey,keysID=keysID, comboFunNwIntNw=comboFunNwIntNw, funNwIntNwID=funNwIntNwID, funNwIntNwText=funNwIntNwText, weightedCheck=weightedCheck, keepAllNodesCheck=keepAllNodesCheck, comboPlotOutput=comboPlotOutput,plotOutputID=plotOutputID, saveCheck=saveCheck,
    plotTitleText=plotTitleText, legendPrefixText=legendPrefixText, showLegendCheck=showLegendCheck, 
    vSizeText=vSizeText, vLabelCexText=vLabelCexText, bgTransparencyText=bgTransparencyText, eColorText=eColorText)
    
    
    # SignalConnects
    RGtk2::gSignalConnect(buttonNetwork, "clicked", generateNetwork, data=list(parentWindow=mainWindow, statusbar=statusbar, 
                                                            commonFields=commonFields, onlyNwFields=fields))
    
    RGtk2::gSignalConnect(comboFunNwIntNw, "changed", selectNetworkType, data=list(comboFunNwIntNw=comboFunNwIntNw,funNwIntNwID=funNwIntNwID, funNwIntNwText=funNwIntNwText, plotTitleText=plotTitleText, comboPlotOutput=comboPlotOutput,plotOutputID=plotOutputID,
                                                                            weightedCheck=weightedCheck, keepAllNodesCheck=keepAllNodesCheck))
    
    return(list(tabPlotNetwork=tabPlotNetwork, fields=fields))
}
