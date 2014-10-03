################################################################################
### Aux functions
# generateReport

# argsList=
# parentWindow, statusbar, 
# expressionText, comboExpressionPlot,
# commonFields$comboFeaTool, commonFields$feaResultsText, 
# commonFields$thresholdAttributeNetworkText, commonFields$thresholdValueNetworkText, commonFields$thresholdOperatorNetworkText
# plotChecksValues

generateReport <- function(button, argsList) # button: Not used, but required for signalConnect.
{ 
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    inputFile <- argsList$commonFields$feaResultsText$getText()
    reportTool <- argsList$commonFields$comboFeaTool$getActive()
    reportTool <- rownames(FEA_tools)[which(as.numeric(FEA_tools[,"ID"])==reportTool)]
    
    # Recargar siempre? donde poner para q solo cuando cambia?
    feaResults <- readFEAresults(feaResultsText=argsList$commonFields$feaResultsText,comboFeaTool=argsList$commonFields$comboFeaTool, expressionText=argsList$commonFields$expressionText,
                                  serverWebReportText=argsList$commonFields$serverWebReportText,alreadyDownloadedCheck=argsList$commonFields$alreadyDownloadedCheck, parentWindow=argsList$parentWindow)
    
    if(is.null(feaResults) || (reportTool==-1 || inputFile==""))
    {
        # Dialog
        dialog <- gtkDialogNewWithButtons("FEA results required", argsList$parentWindow,
                                          c("modal", "destroy-with-parent"), 
                                          "gtk-ok", GtkResponseType["accept"],
                                          show=TRUE)
        dialog[["vbox"]]$add(gtkLabel("Please provide a file with FEA results to continue."))
        response <- dialog$run() 
        dialog$destroy()
        
    }else
    {
        if(grepl(".RData", inputFile, fixed=TRUE)) 
        {
            load(inputFile)
        }
        
        ###################################################################
        # Get arguments
        plotKeggPw <- argsList$plotChecksValues[["plotKeggPw"]]$active 
        plotGoTree <- argsList$plotChecksValues[["plotGoTree"]]$active 
        onlyGoLeaves <- argsList$plotChecksValues[["onlyGoLeaves"]]$active # Not for GtLinker. Only david (&other?)
                                                
        geneExpr <- processExpressionField(argsList$commonFields$expressionText)
        plotExpression <- names(argsList$commonFields$exprPlotOptionsID)[which(argsList$commonFields$exprPlotOptionsID==argsList$commonFields$comboExpressionPlot$getActive())]
        
        filterAttribute <- argsList$commonFields$thresholdAttributeNetworkText$getText()  # Only david and gTLinker
        if(filterAttribute=="") filterAttribute <- NULL
        filterThreshold <- argsList$commonFields$thresholdValueNetworkText$getText()
        if(filterThreshold=="") filterThreshold <- 0
        filterOperator <- argsList$commonFields$thresholdOperatorNetworkText$getText() 
        if(filterOperator=="") filterOperator <- NULL
                    
        # Report
        argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Generating report...")


        tryCatch( 
        {
            FGNet_report(feaResults, geneExpr=geneExpr, plotExpression=plotExpression, onlyGoLeaves=onlyGoLeaves, plotGoTree=plotGoTree, plotKeggPw=plotKeggPw, 
                         filterAttribute=filterAttribute, filterOperator=filterOperator, filterThreshold=filterThreshold)
            argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Report ready.")
        }, error = function (e)
        {
            argsList$statusbar$push(argsList$statusbar$getContextId("info"), "Error generating report. See R console for details.")
            setwd("..")
            stop(e)
        })
        
    }
}

#################################################################################
### Fills Report tab
#################################################################################

tabReport_fill <- function(mainWindow, statusbar, commonFields)
{
    # report_david:
    #   geneLabels=NULL
    
    tabReport <- gtkVBox(FALSE,3)
    
    # plotChecks
    plotChecks <- c("plotKeggPw","plotGoTree", "onlyGoLeaves")  # onlyGoLeaves (only David/other?) TODO
    plotChecksFrame <- gtkFrame("Plot options")
    plotChecksArea <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    plotChecksFrame$add(plotChecksArea)
    plotChecksValues <-list()
    for(ck in plotChecks)
    {
        plotChecksValues[[ck]] <- gtkCheckButton(ck)
        plotChecksValues[[ck]]$active <- TRUE
        plotChecksArea$add(plotChecksValues[[ck]])
    }
    tabReport$packStart(plotChecksFrame, expand=FALSE)
    
    # Report tab end
    buttonReport <- gtkButton("Create")
    tabReport$packStart(buttonReport, expand = FALSE, FALSE, 10) 
    gSignalConnect(buttonReport, "clicked", generateReport, data=list(parentWindow=mainWindow, statusbar=statusbar, 
                                                                      commonFields=commonFields, plotChecksValues=plotChecksValues))               
    return(list(tabReport=tabReport, fields=list(plotChecksValues=plotChecksValues)))
}