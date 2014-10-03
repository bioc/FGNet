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
    tabHelp <- gtkVBox(FALSE,3)
    
    
    
    
    vignetteFrame <- gtkFrame("To see FGNet tutorial, type in R console:")
    gtkFrameSetShadowType(vignetteFrame, GtkShadowType["none"])
    vignetteText <- gtkEntryNew()
    vignetteText$setText('browseVignettes("FGNet")')
    vignetteFrame$add(vignetteText)
    
    tabHelp$packStart(vignetteFrame, FALSE, FALSE, 0)
    
#     buttonVignette <- gtkLabel("User guide")
#     tabHelp$packStart(buttonVignette, FALSE, FALSE, 0)
#     gSignalConnect(buttonVignette, "clicked", showVignette)
# tabHelp$packStart(gtkLinkButtonNewWithLabel("/home/saraa/Desktop/R-refcard.pdf", label = "[ TEEEST ]", show = TRUE), FALSE, FALSE, 0)
        
    #######################################################################
    #  ready
    return(list(tabHelp=tabHelp)) 
    #######################################################################
}
