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
    grTypes <- tolower(FEA_tools[!is.na(argsList$FEA_tools[,"GroupType"]),"GroupType"]) # for gsub
    subGrType <- unique(grTypes[sapply(grTypes, function(x) grepl(x, argsList$thresholdFrame$label))])
    
    argsList$thresholdFrame$label <- gsub(subGrType, tolower(argsList$FEA_tools[reportTool,"GroupType"]), argsList$thresholdFrame$label)
    argsList$thresholdFrame$"tooltip-text" <- gsub(capitalize(subGrType), argsList$FEA_tools[reportTool,"GroupType"], argsList$thresholdFrame$"tooltip-text")
    
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
        
        argsList$frameFEAresultsText$label <- "Job ID or file"        
        
        RGtk2::gtkWidgetSetSensitive(argsList$gtLinkerDownloadExpander, TRUE)
        argsList$plotChecksValues[["plotKeggPw"]]$active <- FALSE
        RGtk2::gtkWidgetSetSensitive(argsList$plotChecksValues[["plotKeggPw"]], FALSE)
        RGtk2::gtkWidgetSetSensitive(argsList$plotChecksValues[["onlyGoLeaves"]], FALSE)
        
    }else{
        argsList$frameFEAresultsText$label <- "FEA results"
        
        argsList$gtLinkerDownloadExpander$Hide()
        RGtk2::gtkWidgetSetSensitive(argsList$gtLinkerDownloadExpander, FALSE)
        RGtk2::gtkWidgetSetSensitive(argsList$plotChecksValues[["plotKeggPw"]], TRUE)
        RGtk2::gtkWidgetSetSensitive(argsList$plotChecksValues[["onlyGoLeaves"]], TRUE)
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
    tabsNw <- RGtk2::gtkNotebookNew()
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
    RGtk2::gtkNotebookInsertPage(tabsNw, tabReport$tabReport, tab.label = gtkLabelNew("HTML report"))
    RGtk2::gtkNotebookInsertPage(tabsNw, tabPlotNetwork$tabPlotNetwork, tab.label = gtkLabelNew("Plot network"))
    # gtkNotebookInsertPage(tabsNw, tabSubNetwork$tabSubNetwork, tab.label = gtkLabelNew("Subnetwork [in construction]"))
    tabNetwork$fields$commonBox$add(tabsNw)
    
    # SignalConnects
    RGtk2::gSignalConnect(tabNetwork$fields$buttonSelectFeaResults, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, feaResultsText=tabNetwork$fields$feaResultsText, comboFeaTool=tabNetwork$fields$comboFeaTool, expressionText=tabNetwork$fields$expressionText, alreadyDownloadedCheck=tabNetwork$fields$alreadyDownloadedCheck)) # Assigns value to feaResultsText_other
    RGtk2::gSignalConnect(tabNetwork$fields$comboFeaTool, "changed", selectNetworkTool, data=c(tabNetwork$fields, list(legendPrefixText=tabPlotNetwork$fields$legendPrefixText, plotChecksValues=tabReport$fields$plotChecksValues, FEA_tools=FEA_tools)))

    fields <- list(commonFields=tabNetwork$fields, tabReport=tabReport, tabPlotNetwork=tabPlotNetwork, tabSubNetwork=tabSubNetwork)
    #######################################################################
    #  ready
    return(c(tabNetwork=tabNetwork$tabNetwork, fields))#, queryArgs=list(comboIdType=comboIdType, annotsArea=annotsArea, emailText=emailText, frameEmail=frameEmail, argsText=argsText)))
    #######################################################################
}
