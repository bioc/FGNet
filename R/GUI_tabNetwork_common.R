
#################################################################################
### Fills the common arguments
#################################################################################

tabNetwork_common_fill <- function(mainWindow,geneList)
{
    tabNetwork <- RGtk2::gtkHBox(FALSE,3)
    
    ############################################################################
    ##### 1. Genes EXPRESSION text field
    # Expression
    table.Left <- RGtk2::gtkTable(rows=20,columns=2)
    frameExpressionFields <- RGtk2::gtkFrame("Expression (only for ploting)")
    RGtk2::gtkFrameSetShadowType(frameExpressionFields, GtkShadowType["none"])
    table.Left$attach(frameExpressionFields, left.attach = 0,2, top.attach = 0,19, xoptions =  c("expand", "fill"), yoptions = c("expand", "fill"))
    tabNetwork$packStart(table.Left, TRUE, TRUE, 10)
    
    # Scrolled area:
    expressionArea <- RGtk2::gtkScrolledWindow()
    expressionArea$setPolicy("automatic", "automatic")
    expressionArea$setShadowType("in")
    frameExpressionFields$add(expressionArea)
    
    # Text field
    expressionText <- RGtk2::gtkTextViewNewWithBuffer()
    expressionText$grabFocus()
    expressionText$"tooltip-text" <- "One gene per line, followed by its expression separated by tab or space"
    expressionText$modifyFont(pangoFontDescriptionFromString(", monospace"))
    expressionArea$add(expressionText)
    
    # Add genes from argument to field
    if(!is.null(geneList) && !is.null(names(geneList)))
    {
        genesBuffer <- RGtk2::gtkTextBufferNew()
        RGtk2::gtkTextBufferSetText(genesBuffer,paste(paste(names(geneList), geneList,sep="\t"),collapse="\n"))
        RGtk2::gtkTextViewSetBuffer(expressionText, genesBuffer)
    }
    
    # Combo
    exprPlotOptionsID <- setNames(0:1, c("border","fill"))
    exprPlotOptionsShow <- setNames(0:1, c("Border","Fill node"))
    frameExpressionPlot <- RGtk2::gtkFrame("Plot as")
    RGtk2::gtkFrameSetShadowType(frameExpressionPlot, GtkShadowType["none"])
    comboExpressionPlot <- RGtk2::gtkComboBoxNewText()
    for (i in capitalize(names(exprPlotOptionsShow))) RGtk2::gtkComboBoxAppendText(comboExpressionPlot, i)
    RGtk2::gtkComboBoxSetTitle(comboExpressionPlot, "comboExpressionPlot")
    comboExpressionPlot$setActive(exprPlotOptionsID["border"])
    comboExpressionPlot$"tooltip-text" <- "Genes with positive values are plotted RED, negative values GREEN"
    frameExpressionPlot$add(comboExpressionPlot)
    table.Left$attach(frameExpressionPlot, left.attach = 0,1, top.attach = 19,20, xoptions = "", yoptions = "")
    
    frameButtonSelectExprFile <- RGtk2::gtkFrameNew("")
    RGtk2::gtkFrameSetShadowType(frameButtonSelectExprFile, GtkShadowType["none"])
    buttonSelectExprFile <- RGtk2::gtkButton("Load from file")   #.RData containing vector or .csv???
    buttonSelectExprFile$"tooltip-text" <- "Text file with one gene per line, followed by its expression separated by TAB"
    buttonSelectExprFile$name <- "buttonSelectExprFile"
    gSignalConnect(buttonSelectExprFile, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, expressionText=expressionText))
    
    frameButtonSelectExprFile$add(buttonSelectExprFile)
    table.Left$attach(frameButtonSelectExprFile, left.attach = 1,2, top.attach = 19,20, xoptions = "", yoptions = "")
        
    ############################################################################
    ##### 2. FGNet options
    
    box.Right <- RGtk2::gtkVBoxNew(FALSE, 10)
    frame.FGNetOptions <- RGtk2::gtkFrameNew("Functional Network")
    #RGtk2::gtkFrameSetShadowType(frame.FGNetOptions, GtkShadowType["none"])
    box.Right$add(frame.FGNetOptions)
    tabNetwork$add(box.Right)
    
    ##### TOP (common box)
    commonBox <- RGtk2::gtkVBoxNew(FALSE, 10)
    commonBox$setBorderWidth(10)
    
    fgnCommon1Box <- RGtk2::gtkHBox(FALSE,3)
    ### Tool
    frameFeaTool <- RGtk2::gtkFrame("FEA tool")
    RGtk2::gtkFrameSetShadowType(frameFeaTool, GtkShadowType["none"])
    frameFeaTool2 <- RGtk2::gtkFrame("")
    RGtk2::gtkFrameSetShadowType(frameFeaTool2, GtkShadowType["none"])
    comboFeaTool <- RGtk2::gtkComboBoxNewText()
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    for (i in rownames(FEA_tools)) RGtk2::gtkComboBoxAppendText(comboFeaTool, i)
    RGtk2::gtkComboBoxSetTitle(comboFeaTool, "comboFeaTool")
    comboFeaTool$setActive(FEA_tools["Imported text file","ID"])
    frameFeaTool$add(comboFeaTool)
    fgnCommon1Box$packStart(frameFeaTool, expand=FALSE)

    # Result
    frameFEAresultsText <- RGtk2::gtkFrame("FEA results")
    RGtk2::gtkFrameSetShadowType(frameFEAresultsText, GtkShadowType["none"])
    feaResultsText <- RGtk2::gtkEntryNew()
    #feaResultsText$setWidthChars(25)
    frameFEAresultsText$add(feaResultsText)
    fgnCommon1Box$packStart(frameFEAresultsText, expand=TRUE)
    
    
#gSignalConnect(feaResultsText, "changed", loadFEAresults?)
    
    # Button
    frameFEAfile <- RGtk2::gtkFrame("")
    RGtk2::gtkFrameSetShadowType(frameFEAfile, GtkShadowType["none"])
    buttonSelectFeaResults <- RGtk2::gtkButton("Select file")
    buttonSelectFeaResults$name <- "buttonSelectFeaResults"
    frameFEAfile$add(buttonSelectFeaResults)
    fgnCommon1Box$packStart(frameFEAfile, expand=FALSE)
    
    commonBox$packStart(fgnCommon1Box, expand=FALSE)
                    
    ### GeneTerm Linker:
    gtLinkerDownloadExpander <- RGtk2::gtkExpander("GeneTerm Linker download options", show=FALSE)
    gtLinkerDownloadExpander$Hide()
    gtLinkerDownloadBox <- RGtk2::gtkHBox(FALSE,3)
    
    # Server
    serverWebReportFrame <- RGtk2::gtkFrame("Server")
    RGtk2::gtkFrameSetShadowType(serverWebReportFrame, GtkShadowType["none"])
    serverWebReportText <- RGtk2::gtkEntryNew()
    serverWebReportText$setWidthChars(25)
    serverWebReportText$setText("http://gtlinker.cnb.csic.es")
    serverWebReportFrame$add(serverWebReportText)
    gtLinkerDownloadBox$packStart(serverWebReportFrame, expand=FALSE)
    # Checkbox
    alreadyDownloadedCheckFrame <- RGtk2::gtkFrame("")
    RGtk2::gtkFrameSetShadowType(alreadyDownloadedCheckFrame, GtkShadowType["none"])
    alreadyDownloadedCheck<- RGtk2::gtkCheckButton("AlreadyDownloaded")
    alreadyDownloadedCheckFrame$add(alreadyDownloadedCheck)
    gtLinkerDownloadBox$packStart(alreadyDownloadedCheckFrame, expand=FALSE)
    # Add
    gtLinkerDownloadExpander$add(gtLinkerDownloadBox)
    RGtk2::gtkWidgetSetSensitive(gtLinkerDownloadExpander, FALSE)
    gtLinkerDownloadExpander$Hide()
    commonBox$packStart(gtLinkerDownloadExpander, expand=FALSE)
    
    fgnCommon2Box <- RGtk2::gtkHBox(FALSE,3)
    
    
    ### Threshold frame
    thresholdFrame <- RGtk2::gtkFrame("Filter clusters") # MG: silhouette, CL... # OTHER:  TODO (desactivar si no es david/gtlinker?)
    thresholdFrame$"tooltip-text" <- "Clusters that meet these criteria will be REMOVED from the network"
    thresholdBoxNetwork <- RGtk2::gtkHBox(FALSE,3)
    
    thresholdAttributeNetworkFrame <- RGtk2::gtkFrame("Attribute")
    RGtk2::gtkFrameSetShadowType(thresholdAttributeNetworkFrame, GtkShadowType["none"])
    thresholdAttributeNetworkText <- RGtk2::gtkEntryNew()
    thresholdAttributeNetworkText$setWidthChars(15)
    thresholdAttributeNetworkText$setText("")
    thresholdAttributeNetworkFrame$add(thresholdAttributeNetworkText)
    thresholdBoxNetwork$packStart(thresholdAttributeNetworkFrame, expand=TRUE,10)
    
    thresholdOperatorNetworkFrame <- RGtk2::gtkFrame("Operator") 
    RGtk2::gtkFrameSetShadowType(thresholdOperatorNetworkFrame, GtkShadowType["none"])
    thresholdOperatorNetworkText <- RGtk2::gtkEntryNew()
    thresholdOperatorNetworkText$setWidthChars(5)
    thresholdOperatorNetworkText$setText("<")
    thresholdOperatorNetworkFrame$add(thresholdOperatorNetworkText)
    thresholdBoxNetwork$packStart(thresholdOperatorNetworkFrame, expand=FALSE,10)
    
    thresholdValueNetworkFrame <- RGtk2::gtkFrame("Value") 
    RGtk2::gtkFrameSetShadowType(thresholdValueNetworkFrame, GtkShadowType["none"])
    thresholdValueNetworkText <- RGtk2::gtkEntryNew()
    thresholdValueNetworkText$setWidthChars(5)
    thresholdValueNetworkText$setText("0")
    thresholdValueNetworkFrame$add(thresholdValueNetworkText)
    thresholdBoxNetwork$packStart(thresholdValueNetworkFrame, expand=FALSE,10)
    thresholdFrame$add(thresholdBoxNetwork)
    fgnCommon2Box$packStart(thresholdFrame, expand=TRUE)
    commonBox$packStart(fgnCommon2Box, expand=FALSE)

    frame.FGNetOptions$add(commonBox)

    
    # List fields:
    fields <- list(expressionText=expressionText, comboExpressionPlot=comboExpressionPlot, exprPlotOptionsID=exprPlotOptionsID, buttonSelectExprFile=buttonSelectExprFile, 
                   commonBox=commonBox,
                   comboFeaTool=comboFeaTool, frameFEAresultsText=frameFEAresultsText, feaResultsText=feaResultsText,
                   serverWebReportText=serverWebReportText,alreadyDownloadedCheck=alreadyDownloadedCheck,
                   thresholdFrame=thresholdFrame, thresholdAttributeNetworkText=thresholdAttributeNetworkText, 
                   thresholdOperatorNetworkText=thresholdOperatorNetworkText, thresholdValueNetworkText=thresholdValueNetworkText, 
                   gtLinkerDownloadExpander=gtLinkerDownloadExpander, buttonSelectFeaResults=buttonSelectFeaResults)

    return(list(tabNetwork=tabNetwork, fields=fields))
}
