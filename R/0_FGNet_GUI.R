FGNet_GUI <- function(geneList=NULL)
{
    if(!loadInstPkg("RGtk2")) stop("Package RGtk2 is required to use the GUI.")
    
    # feaResults <<- NULL # Global

    data("FEA_tools", envir = environment())
    FEA_tools <- get("FEA_tools", envir  = environment())
    
    data("groupTypes", envir = environment())
    groupTypes <- get("groupTypes", envir  = environment())
    
    data("organisms", envir = environment())
    organisms <- get("organisms", envir  = environment())
    
    data("GOEvidenceCodes", envir = environment())
    GOEvidenceCodes <- get("GOEvidenceCodes", envir  = environment())
    
    # Crear ventana
    mainWindow <- gtkWindowNew("toplevel", show=FALSE)
    #mainWindow$setDefaultSize(500,768)
    mainWindow$setTitle("Functional Gene Networks")
    mainWindow$setBorderWidth(2)
    
    # Statusbar
    statusbar <- gtkStatusbar()
    statusbar$push(statusbar$getContextId("info"), "Ready.")
    
    # Create tabs
    tabsSteps <- gtkNotebookNew()
    tabsSteps$name<-"tabsSteps"
    tabsSteps$setTabPos("left")
    
    
    ###################################################################################
    ############## Network tab
    tabNetwork <- tabNetwork_fill(mainWindow, statusbar, geneList)
    
    ###################################################################################
    ############## FEA tab
    tabFEA <- tabFEA_fill(mainWindow, statusbar, geneList)
    
    ###################################################################################
    ############## Other tabs
    # (tabCompare)
    # TO DO
    
    # (tabHelp)
    tabHelp <- tabHelp_fill(mainWindow) # TO DO
    
    ###################################################################################
    ############## Add tabs to mainWindow
    tabID <- setNames(0:4, c("FEA", "Network", "Compare", "help")) #"Import",
    gtkNotebookInsertPage(tabsSteps, tabFEA$tabFEA, position = tabID["FEA"], tab.label = gtkLabelNew("1 - FEA    "))
    # gtkNotebookInsertPage(tabsSteps, tabImport, position = tabID["Import"], tab.label = gtkLabelNew("import?"))
    gtkNotebookInsertPage(tabsSteps, tabNetwork$tabNetwork, position = tabID["Network"], tab.label =  gtkLabelNew("2 - Network"))
    # gtkNotebookInsertPage(tabsSteps, tabCompare, position = tabID["Compare"], tab.label = gtkLabelNew("compare?"))
    gtkNotebookInsertPage(tabsSteps, tabHelp$tabHelp, position = tabID["help"], tab.label = gtkLabelNew("Help"))

    box.vMain <- gtkVBoxNew()
    box.vMain$packStart(tabsSteps, expand=TRUE, fill=TRUE, 0)
    box.vMain$packStart(statusbar, expand=FALSE, fill=TRUE, 0)
    

    ###################################################################################
    ############## Signal connects
    # Submit Query:
    gSignalConnect(tabFEA$commonFields$buttonSubmit, "clicked", submitQuery, data=list(parentWindow=mainWindow, statusbar=statusbar, 
                                                                                       feaResultsText=tabNetwork$commonFields$feaResultsText, comboFeaTool=tabNetwork$commonFields$comboFeaTool, expressionText=tabNetwork$commonFields$expressionText,serverWebReportText=tabNetwork$commonFields$serverWebReportText, 
                                                                                       tabsSteps=tabsSteps, tabID=tabID, tabFEA=tabFEA))    
    mainWindow$add(box.vMain)
    mainWindow$show()
}
