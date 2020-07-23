################################################################################
### Aux functions

showVignette <- function(button) # button: Not used, but required for signalConnect.
{
    #vignette("UserGuide", "FGNet")
    #openPDF(file, bg=TRUE)
    #print('vignette("builtinMethods", "affy")')
}

#################################################################################
### Fills 
#################################################################################
 
tabHelp_fill <- function(mainWindow)
{
    tabHelp <- RGtk2::gtkHBox(FALSE,10)
    tabHelpV <- RGtk2::gtkVBox(FALSE,10)
    
    helpTxt <- 'To generate the functional network:\n\n 1. Run the Functional Enrichment Analisis (FEA) with one of the tools (Tab "1 - FEA")\n\n 2. Generate the HTML report or personalize the networks (Tab "2 - Network")\n'
    tabHelpV$packStart(RGtk2::gtkLabel(helpTxt), TRUE, TRUE, 50)
    
    vignetteFrame <- RGtk2::gtkFrame("To see FGNet tutorial, type in R console:")
    RGtk2::gtkFrameSetShadowType(vignetteFrame, GtkShadowType["none"])
    vignetteText <- RGtk2::gtkEntryNew()
    vignetteText$setText('browseVignettes("FGNet")')
    vignetteFrame$add(vignetteText)
    
    tabHelpV$packStart(vignetteFrame, FALSE, FALSE, 50)
    
    tabHelp$packStart(tabHelpV, FALSE, FALSE, 50)
    
#     buttonVignette <- RGtk2::gtkLabel("User guide")
#     tabHelp$packStart(buttonVignette, FALSE, FALSE, 0)
#     gSignalConnect(buttonVignette, "clicked", showVignette)
# tabHelp$packStart(RGtk2::gtkLinkButtonNewWithLabel("/home/saraa/Desktop/R-refcard.pdf", label = "[ TEEEST ]", show = TRUE), FALSE, FALSE, 0)
        
    #######################################################################
    #  ready
    return(list(tabHelp=tabHelp)) 
    #######################################################################
}
