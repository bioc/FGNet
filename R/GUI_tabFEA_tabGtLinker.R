#################################################################################
### Fills tabGTL
#################################################################################

tabGTL_fill <- function(mainWindow)
{
    ########################################################################
    # Fill tab
    tabGTL <- RGtk2::gtkVBox(homogeneous=FALSE, spacing=3)
        
    # Organism
    comboOrg <- RGtk2::gtkComboBoxNewText()
    GtL_org <- as.character(0:1)
    names(GtL_org) <- c("Hs","Sc")
    GtL_org_full <- c("Homo sapiens","Saccharomyces cerevisiae")
    names(GtL_org_full) <- c("Hs","Sc")
    
    for (ch in GtL_org_full) RGtk2::gtkComboBoxAppendText(comboOrg, ch)
    comboOrg$setActive(GtL_org["Hs"])
    RGtk2::gtkComboBoxSetTitle(comboOrg, "Org")
    frameOrg <- RGtk2::gtkFrame("Organism")
    RGtk2::gtkFrameSetShadowType(frameOrg, GtkShadowType["none"])
    frameOrg$add(comboOrg)
    tabGTL$packStart(frameOrg, expand=FALSE)
    
    # Annotations
    GtL_annots <- c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "InterPro_Motifs")
    frameAnnots <- RGtk2::gtkFrame("Annotations")
    annotsArea <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=0)
    frameAnnots$add(annotsArea)
    checkAnnotsGtL<-list()
    for(annot in GtL_annots)
    {
      checkAnnotsGtL[[annot]] <- RGtk2::gtkCheckButton(annot)
      checkAnnotsGtL[[annot]]$active <- TRUE
      annotsArea$add(checkAnnotsGtL[[annot]])
    }
    tabGTL$packStart(frameAnnots, expand=FALSE)
    
    # MinSupport
    comboMinSupport <- RGtk2::gtkComboBoxNewText()
    GtL_minSupport <- as.character(0:7)
    names(GtL_minSupport) <- as.character(2:9)
    for (ch in names(GtL_minSupport)) RGtk2::gtkComboBoxAppendText(comboMinSupport, ch) 
    comboMinSupport$setActive(GtL_minSupport["4"])
    comboMinSupport$"tooltip-text" <- "Minimum size of the gene-term set\n(in number of genes)"
    RGtk2::gtkComboBoxSetTitle(comboMinSupport, "MinSuport")
    frameMinSupport <- RGtk2::gtkFrame("Minimum support")
    RGtk2::gtkFrameSetShadowType(frameMinSupport, GtkShadowType["none"])
    frameMinSupport$add(comboMinSupport)
    tabGTL$packStart(frameMinSupport, expand=FALSE)
    #tabGTL$packStart(comboMinSupport)
    
    # serverWeb
    frameServerWeb <- RGtk2::gtkFrame("Server")
    RGtk2::gtkFrameSetShadowType(frameServerWeb, GtkShadowType["none"])
    serverWebText <- RGtk2::gtkEntryNew()
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
