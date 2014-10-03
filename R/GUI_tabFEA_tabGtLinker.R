#################################################################################
### Fills tabGTL
#################################################################################

tabGTL_fill <- function(mainWindow)
{
    ########################################################################
    # Fill tab
    tabGTL <- gtkVBox(homogeneous=FALSE, spacing=3)
        
    # Organism
    comboOrg <- gtkComboBoxNewText()
    GtL_org <- as.character(0:1)
    names(GtL_org) <- c("Hs","Sc")
    GtL_org_full <- c("Homo sapiens","Saccharomyces cerevisiae")
    names(GtL_org_full) <- c("Hs","Sc")
    
    for (ch in GtL_org_full) gtkComboBoxAppendText(comboOrg, ch)
    comboOrg$setActive(GtL_org["Hs"])
    gtkComboBoxSetTitle(comboOrg, "Org")
    frameOrg <- gtkFrame("Organism")
    gtkFrameSetShadowType(frameOrg, GtkShadowType["none"])
    frameOrg$add(comboOrg)
    tabGTL$packStart(frameOrg, expand=FALSE)
    
    # Annotations
    GtL_annots <- c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs")
    frameAnnots <- gtkFrame("Annotations")
    annotsArea <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    frameAnnots$add(annotsArea)
    checkAnnotsGtL<-list()
    for(annot in GtL_annots)
    {
      checkAnnotsGtL[[annot]] <- gtkCheckButton(annot)
      checkAnnotsGtL[[annot]]$active <- TRUE
      annotsArea$add(checkAnnotsGtL[[annot]])
    }
    tabGTL$packStart(frameAnnots, expand=FALSE)
    
    # MinSupport
    comboMinSupport <- gtkComboBoxNewText()
    GtL_minSupport <- as.character(0:7)
    names(GtL_minSupport) <- as.character(2:9)
    for (ch in names(GtL_minSupport)) gtkComboBoxAppendText(comboMinSupport, ch) 
    comboMinSupport$setActive(GtL_minSupport["4"])
    gtkComboBoxSetTitle(comboMinSupport, "MinSuport")
    frameMinSupport <- gtkFrame("Minimum support")
    gtkFrameSetShadowType(frameMinSupport, GtkShadowType["none"])
    frameMinSupport$add(comboMinSupport)
    tabGTL$packStart(frameMinSupport, expand=FALSE)
    #tabGTL$packStart(comboMinSupport)
    
    # serverWeb
    frameServerWeb <- gtkFrame("Server")
    gtkFrameSetShadowType(frameServerWeb, GtkShadowType["none"])
    serverWebText <- gtkEntryNew()
    serverWebText$setWidthChars(25)
    serverWebText$setText("http://gtlinker.cnb.csic.es")
    frameServerWeb$add(serverWebText)
    tabGTL$packStart(frameServerWeb, expand=FALSE)


    fields <- list(comboOrg=comboOrg, GtL_org=GtL_org,
               checkAnnotsGtL=checkAnnotsGtL, GtL_annots=GtL_annots,
               comboMinSupport=comboMinSupport, GtL_minSupport=GtL_minSupport, serverWebText=serverWebText)

    #######################################################################
    ## tabGTL ready    
    return(list(tabGTL=tabGTL, queryArgs=fields))
    #######################################################################
}