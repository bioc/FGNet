################################################################################
### Aux functions
# selectNetworkTool
selectNetworkTool <- function(aux, argsList) 
{
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    data("groupTypes", envir = environment())
    groupTypes<- get("groupTypes", envir  = environment())
    
    reportTool <- rownames(argsList$FEA_tools[which(as.numeric(argsList$FEA_tools[,"ID"])==argsList$comboFeaTool$getActive()),, drop=FALSE])
    grTypes <- FEA_tools[!is.na(argsList$FEA_tools[,"GroupType"]),"GroupType"] # for gsub
    subGrType <- unique(grTypes[sapply(grTypes, function(x) grepl(x, argsList$thresholdFrame$label))])
    
    argsList$thresholdFrame$label <- gsub(subGrType, argsList$FEA_tools[reportTool,"GroupType"], argsList$thresholdFrame$label)
    
    argsList$thresholdAttributeNetworkText$setText(argsList$FEA_tools[reportTool,"DefaultFilter"])
    if(argsList$thresholdAttributeNetworkText$getText()=="NA") argsList$thresholdAttributeNetworkText$setText("")
    argsList$thresholdOperatorNetworkText$setText(argsList$FEA_tools[reportTool,"DefaultFiltOperator"])
    if(argsList$thresholdOperatorNetworkText$getText()=="NA") argsList$thresholdOperatorNetworkText$setText("")
    argsList$thresholdValueNetworkText$setText(argsList$FEA_tools[reportTool,"DefaultFiltThreshold"])
    if(argsList$thresholdValueNetworkText$getText()=="NA") argsList$thresholdValueNetworkText$setText("")
    #gtkWidgetSetSensitive(thresholdFrame, !is.na(argsList$FEA_tools[reportTool,"DefaultFilter"]))
    
    argsList$legendPrefixText$setText(groupTypes[argsList$FEA_tools[reportTool,"GroupType"], "prefix"])        
    
    if(reportTool == "GeneTerm Linker")
    {        
        argsList$gtLinkerDownloadExpander$Show()
        gtkWidgetSetSensitive(argsList$gtLinkerDownloadExpander, TRUE)
        argsList$plotChecksValues[["plotKeggPw"]]$active <- FALSE
        gtkWidgetSetSensitive(argsList$plotChecksValues[["plotKeggPw"]], FALSE)
        gtkWidgetSetSensitive(argsList$plotChecksValues[["onlyGoLeaves"]], FALSE)
        
    }else{
        argsList$gtLinkerDownloadExpander$Hide()
        gtkWidgetSetSensitive(argsList$gtLinkerDownloadExpander, FALSE)
        gtkWidgetSetSensitive(argsList$plotChecksValues[["plotKeggPw"]], TRUE)
        gtkWidgetSetSensitive(argsList$plotChecksValues[["onlyGoLeaves"]], TRUE)
    }    
}

#################################################################################
### Fills network tab
#################################################################################

tabNetwork_fill <- function(mainWindow, statusbar, geneList) 
{
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    ####################################################
    ####### Common
    tabNetwork <- tabNetwork_common_fill(mainWindow, geneList)

    # Create subtabs: 
    tabsNw <- gtkNotebookNew()
    tabsNw$setTabPos("top")
   
    ####################################################
    ####### Report tab

    tabReport <- tabReport_fill(mainWindow, statusbar, commonFields=tabNetwork$fields)
 
    ####################################################
    ####### Only network tab
    
    # clustersId  # Global for subnetwork
    tabPlotNetwork <- tabPlotNetwork_fill(mainWindow, statusbar, commonFields=tabNetwork$fields)

    
    # ####################################################
    # ####### SubNw network tab

    tabSubNetwork <- tabSubNetwork_fill()
    
    # ####################################################
    # ####### FGNet tab end
    gtkNotebookInsertPage(tabsNw, tabReport$tabReport, tab.label = gtkLabelNew("HTML report"))
    gtkNotebookInsertPage(tabsNw, tabPlotNetwork$tabPlotNetwork, tab.label = gtkLabelNew("Plot network"))
    # gtkNotebookInsertPage(tabsNw, tabSubNetwork$tabSubNetwork, tab.label = gtkLabelNew("Subnetwork [in construction]"))
    tabNetwork$fields$commonBox$add(tabsNw)
    
    # SignalConnects
    gSignalConnect(tabNetwork$fields$buttonSelectFeaResults, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, feaResultsText=tabNetwork$fields$feaResultsText, comboFeaTool=tabNetwork$fields$comboFeaTool, expressionText=tabNetwork$fields$expressionText, alreadyDownloadedCheck=tabNetwork$fields$alreadyDownloadedCheck)) # Assigns value to feaResultsText_other
    gSignalConnect(tabNetwork$fields$comboFeaTool, "changed", selectNetworkTool, data=c(tabNetwork$fields, list(legendPrefixText=tabPlotNetwork$fields$legendPrefixText, plotChecksValues=tabReport$fields$plotChecksValues, FEA_tools=FEA_tools)))

    fields <- list(commonFields=tabNetwork$fields, tabReport=tabReport, tabPlotNetwork=tabPlotNetwork, tabSubNetwork=tabSubNetwork)
    #######################################################################
    #  ready
    return(c(tabNetwork=tabNetwork$tabNetwork, fields))#, queryArgs=list(comboIdType=comboIdType, annotsArea=annotsArea, emailText=emailText, frameEmail=frameEmail, argsText=argsText)))
    #######################################################################
}
