
#################################################################################
### Fills the common arguments
#################################################################################

tabNetwork_common_fill <- function(mainWindow,geneList)
{
    tabNetwork <- gtkHBox(FALSE,3)
    
    ############################################################################
    ##### 1. Genes EXPRESSION text field
    # Expression
    table.Left <- gtkTable(rows=20,columns=2)
    frameExpressionFields <- gtkFrame("Expression (only for ploting)")
    gtkFrameSetShadowType(frameExpressionFields, GtkShadowType["none"])
    table.Left$attach(frameExpressionFields, left.attach = 0,2, top.attach = 0,19, xoptions =  c("expand", "fill"), yoptions = c("expand", "fill"))
    tabNetwork$packStart(table.Left, TRUE, TRUE, 10)
    
    # Scrolled area:
    expressionArea <- gtkScrolledWindow()
    expressionArea$setPolicy("automatic", "automatic")
    expressionArea$setShadowType("in")
    frameExpressionFields$add(expressionArea)
    
    # Text field
    expressionText <- gtkTextViewNewWithBuffer()
    expressionText$grabFocus()
    expressionText$"tooltip-text" <- "One gene per line, followed by its expression separated by tab or space"
    expressionText$modifyFont(pangoFontDescriptionFromString(", monospace"))
    expressionArea$add(expressionText)
    
    # Add genes from argument to field
    if(!is.null(geneList) && !is.null(names(geneList)))
    {
        genesBuffer <- gtkTextBufferNew()
        gtkTextBufferSetText(genesBuffer,paste(paste(names(geneList), geneList,sep="\t"),collapse="\n"))
        gtkTextViewSetBuffer(expressionText, genesBuffer)
    }
    
    # Combo
    exprPlotOptionsID <- setNames(0:1, c("border","fill"))
    exprPlotOptionsShow <- setNames(0:1, c("Border","Fill node"))
    frameExpressionPlot <- gtkFrame("Plot as")
    gtkFrameSetShadowType(frameExpressionPlot, GtkShadowType["none"])
    comboExpressionPlot <- gtkComboBoxNewText()
    for (i in capitalize(names(exprPlotOptionsShow))) gtkComboBoxAppendText(comboExpressionPlot, i)
    gtkComboBoxSetTitle(comboExpressionPlot, "comboExpressionPlot")
    comboExpressionPlot$setActive(exprPlotOptionsID["border"])
    comboExpressionPlot$"tooltip-text" <- "Genes with positive values are plotted RED, negative values GREEN"
    frameExpressionPlot$add(comboExpressionPlot)
    table.Left$attach(frameExpressionPlot, left.attach = 0,1, top.attach = 19,20, xoptions = "", yoptions = "")
    
    frameButtonSelectExprFile <- gtkFrameNew("")
    gtkFrameSetShadowType(frameButtonSelectExprFile, GtkShadowType["none"])
    buttonSelectExprFile <- gtkButton("Load from file")   #.RData containing vector or .csv???
    buttonSelectExprFile$"tooltip-text" <- "Text file with one gene per line, followed by its expression separated by TAB"
    buttonSelectExprFile$name <- "buttonSelectExprFile"
    gSignalConnect(buttonSelectExprFile, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, expressionText=expressionText))
    
    frameButtonSelectExprFile$add(buttonSelectExprFile)
    table.Left$attach(frameButtonSelectExprFile, left.attach = 1,2, top.attach = 19,20, xoptions = "", yoptions = "")
        
    ############################################################################
    ##### 2. FGNet options
    
    box.Right <- gtkVBoxNew(FALSE, 10)
    frame.FGNetOptions <- gtkFrameNew("Functional Network")
    #gtkFrameSetShadowType(frame.FGNetOptions, GtkShadowType["none"])
    box.Right$add(frame.FGNetOptions)
    tabNetwork$add(box.Right)
    
    ##### TOP (common box)
    commonBox <- gtkVBoxNew(FALSE, 10)
    commonBox$setBorderWidth(10)
    
    fgnCommon1Box <- gtkHBox(FALSE,3)
    ### Tool
    frameFeaTool <- gtkFrame("FEA tool")
    gtkFrameSetShadowType(frameFeaTool, GtkShadowType["none"])
    frameFeaTool2 <- gtkFrame("")
    gtkFrameSetShadowType(frameFeaTool2, GtkShadowType["none"])
    comboFeaTool <- gtkComboBoxNewText()
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    for (i in rownames(FEA_tools)) gtkComboBoxAppendText(comboFeaTool, i)
    gtkComboBoxSetTitle(comboFeaTool, "comboFeaTool")
    comboFeaTool$setActive(FEA_tools["Imported text file","ID"])
    frameFeaTool$add(comboFeaTool)
    fgnCommon1Box$packStart(frameFeaTool, expand=FALSE)

    # Result
    frameFEAresultsText <- gtkFrame("FEA results")
    gtkFrameSetShadowType(frameFEAresultsText, GtkShadowType["none"])
    feaResultsText <- gtkEntryNew()
    #feaResultsText$setWidthChars(25)
    frameFEAresultsText$add(feaResultsText)
    fgnCommon1Box$packStart(frameFEAresultsText, expand=TRUE)
    
    
#gSignalConnect(feaResultsText, "changed", loadFEAresults?)
    
    # Button
    frameFEAfile <- gtkFrame("")
    gtkFrameSetShadowType(frameFEAfile, GtkShadowType["none"])
    buttonSelectFeaResults <- gtkButton("Select file")
    buttonSelectFeaResults$name <- "buttonSelectFeaResults"
    frameFEAfile$add(buttonSelectFeaResults)
    fgnCommon1Box$packStart(frameFEAfile, expand=FALSE)
    
    commonBox$packStart(fgnCommon1Box, expand=FALSE)
                    
    ### GeneTerm Linker:
    gtLinkerDownloadExpander <- gtkExpander("GeneTerm Linker download options", show=FALSE)
    gtLinkerDownloadExpander$Hide()
    gtLinkerDownloadBox <- gtkHBox(FALSE,3)
    
    # Server
    serverWebReportFrame <- gtkFrame("Server")
    gtkFrameSetShadowType(serverWebReportFrame, GtkShadowType["none"])
    serverWebReportText <- gtkEntryNew()
    serverWebReportText$setWidthChars(25)
    serverWebReportText$setText("http://gtlinker.cnb.csic.es")
    serverWebReportFrame$add(serverWebReportText)
    gtLinkerDownloadBox$packStart(serverWebReportFrame, expand=FALSE)
    # Checkbox
    alreadyDownloadedCheckFrame <- gtkFrame("")
    gtkFrameSetShadowType(alreadyDownloadedCheckFrame, GtkShadowType["none"])
    alreadyDownloadedCheck<- gtkCheckButton("AlreadyDownloaded")
    alreadyDownloadedCheckFrame$add(alreadyDownloadedCheck)
    gtLinkerDownloadBox$packStart(alreadyDownloadedCheckFrame, expand=FALSE)
    # Add
    gtLinkerDownloadExpander$add(gtLinkerDownloadBox)
    gtkWidgetSetSensitive(gtLinkerDownloadExpander, FALSE)
    gtLinkerDownloadExpander$Hide()
    commonBox$packStart(gtLinkerDownloadExpander, expand=FALSE)
    
    fgnCommon2Box <- gtkHBox(FALSE,3)
    
    
    ### Threshold frame
    thresholdFrame <- gtkFrame("Filter clusters") # MG: silhouette, CL... # OTHER:  TODO (desactivar si no es david/gtlinker?)
    thresholdFrame$"tooltip-text" <- "Clusters that meet these criteria will be REMOVED from the network"
    thresholdBoxNetwork <- gtkHBox(FALSE,3)
    
    thresholdAttributeNetworkFrame <- gtkFrame("Attribute")
    gtkFrameSetShadowType(thresholdAttributeNetworkFrame, GtkShadowType["none"])
    thresholdAttributeNetworkText <- gtkEntryNew()
    thresholdAttributeNetworkText$setWidthChars(15)
    thresholdAttributeNetworkText$setText("")
    thresholdAttributeNetworkFrame$add(thresholdAttributeNetworkText)
    thresholdBoxNetwork$packStart(thresholdAttributeNetworkFrame, expand=TRUE,10)
    
    thresholdOperatorNetworkFrame <- gtkFrame("Operator") 
    gtkFrameSetShadowType(thresholdOperatorNetworkFrame, GtkShadowType["none"])
    thresholdOperatorNetworkText <- gtkEntryNew()
    thresholdOperatorNetworkText$setWidthChars(5)
    thresholdOperatorNetworkText$setText("<")
    thresholdOperatorNetworkFrame$add(thresholdOperatorNetworkText)
    thresholdBoxNetwork$packStart(thresholdOperatorNetworkFrame, expand=FALSE,10)
    
    thresholdValueNetworkFrame <- gtkFrame("Value") 
    gtkFrameSetShadowType(thresholdValueNetworkFrame, GtkShadowType["none"])
    thresholdValueNetworkText <- gtkEntryNew()
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