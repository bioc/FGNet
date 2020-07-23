################################################################################
### Aux functions

#################################################################################
### Fills the Cluster (subnetwork) tab
#################################################################################


tabSubNetwork_fill <- function()
{
    tabSubNetwork <- RGtk2::gtkVBox(FALSE,3)
    
    
    fgnSubNwVBox <- RGtk2::gtkVBox(FALSE,3)
    fgnSubNwL1Box <- RGtk2::gtkHBox(FALSE,3)
    fgnSubNwL2Box <- RGtk2::gtkHBox(FALSE,3)
    
    # Select cluster
    clSelectClusterFrame <- RGtk2::gtkFrame("Select cluster")
    RGtk2::gtkFrameSetShadowType(clSelectClusterFrame, GtkShadowType["none"])
    clSelectCluster <- RGtk2::gtkHBox(FALSE,3)
    txtSelectCluster <- RGtk2::gtkEntryNew()
    txtSelectCluster$setWidthChars(5)
    comboSelectCluster <- RGtk2::gtkComboBoxNewText(show=FALSE)
    clSelectCluster$packStart(txtSelectCluster, expand=FALSE)
    clSelectCluster$packStart(comboSelectCluster, expand=FALSE)
    clSelectClusterFrame$add(clSelectCluster)
    fgnSubNwL1Box$packStart(clSelectClusterFrame, expand=FALSE)
    
    # plotOutput="static"
    clplotOutputFrame <- RGtk2::gtkFrame("Output")
    RGtk2::gtkFrameSetShadowType(clplotOutputFrame, GtkShadowType["none"])
    comboClplotOutput <- RGtk2::gtkComboBoxNewText()
#for (i in capitalize(names(plotOutputID)))RGtk2::gtkComboBoxAppendText(comboClplotOutput, i)
#comboClplotOutput$setActive(plotOutputID["static"])
    clplotOutputFrame$add(comboClplotOutput)
    fgnSubNwL1Box$packStart(clplotOutputFrame, expand=FALSE)
    
    # Checkbox
    clPlotSaveCheck<- RGtk2::gtkCheckButton("Save (.RData)")
    clPlotSaveCheck$"tooltip-text" <- "Save igraph & matrices"
    fgnSubNwL1Box$packStart(clPlotSaveCheck, expand=FALSE)
    
    fgnSubNwFrame <- RGtk2::gtkFrameNew("Plot...")
    nwVBox1 <- RGtk2::gtkVBox()
    nwButton1 <- RGtk2::gtkButton("nw1")
    nwVBox1$packEnd(nwButton1, expand=FALSE)
#gSignalConnect(nwButton1, "clicked", generateNetwork)
    
    nwVBox2 <- RGtk2::gtkVBox()
    nwButton2 <- RGtk2::gtkButton("nw2")
    nwVBox2$packEnd(nwButton2, expand=FALSE)
#gSignalConnect(nwButton2, "clicked", generateNetwork)
    
    nwVBox3 <- RGtk2::gtkVBox()
    nwButton3 <- RGtk2::gtkButton("nw3")
    nwVBox3$packEnd(nwButton3, expand=FALSE)
# gSignalConnect(nwButton3, "clicked", generateNetwork)
    
    nwVBox5 <- RGtk2::gtkVBox()
    nwButton5 <- RGtk2::gtkButton("nw5")
    nwVBox5$packEnd(nwButton5, expand=FALSE)
#    gSignalConnect(nwButton5, "clicked", generateNetwork)
    
    nwVBox6 <- RGtk2::gtkVBox()
    nwButton6 <- RGtk2::gtkButton("nw6")
    nwVBox6$packEnd(nwButton6, expand=FALSE)
#    gSignalConnect(nwButton6, "clicked", generateNetwork)
    
    nwVBox7 <- RGtk2::gtkVBox()
    nwButton7 <- RGtk2::gtkButton("nw7")
    nwVBox7$packEnd(nwButton7, expand=FALSE)
#    gSignalConnect(nwButton7, "clicked", generateNetwork)
    
    nwVBox8 <- RGtk2::gtkVBox()
    nwButton8 <- RGtk2::gtkButton("nw8")
    nwVBox8$packEnd(nwButton8, expand=FALSE)
#   gSignalConnect(nwButton8, "clicked", generateNetwork)
    
    
    fgnSubNwTable <- RGtk2::gtkTable(rows = 6, columns = 3, homogeneous = FALSE)
    fgnSubNwTable$attach(RGtk2::gtkLabel("(background)"), left.attach=1,2, top.attach=0,1)
    fgnSubNwTable$attach(RGtk2::gtkLabel("(bipartite)"), left.attach=2,3, top.attach=0,1)
    fgnSubNwTable$attach(RGtk2::gtkLabel("Gene-term sets"), left.attach=0,1, top.attach=1,2)
    fgnSubNwTable$attach(RGtk2::gtkLabel(""), left.attach=0,1, top.attach=3,4)
    fgnSubNwTable$attach(RGtk2::gtkLabel("Terms sharing genes"), left.attach=0,1, top.attach=4,5)
    fgnSubNwTable$attach(RGtk2::gtkLabel("Genes sharing terms"), left.attach=0,1, top.attach=5,6)
    
    fgnSubNwTable$attach(nwVBox1, left.attach=1,2, top.attach=1,2)
    fgnSubNwTable$attach(nwVBox2, left.attach=2,3, top.attach=1,2)
    
    fgnSubNwTable$attach(nwVBox3, left.attach=1,2, top.attach=2,3)
    
    fgnSubNwTable$attach(nwVBox5, left.attach=1,2, top.attach=4,5) 
    fgnSubNwTable$attach(nwVBox6, left.attach=2,3, top.attach=4,5) 
    
    fgnSubNwTable$attach(nwVBox7, left.attach=1,2, top.attach=5,6) 
    fgnSubNwTable$attach(nwVBox8, left.attach=2,3, top.attach=5,6) 
    
    fgnSubNwTable$setColSpacing(0, 5)
    
    fgnSubNwFrame$add(fgnSubNwTable)
    fgnSubNwVBox$packStart(fgnSubNwL1Box, expand=FALSE,10)
    fgnSubNwVBox$packStart(fgnSubNwL2Box, expand=FALSE,10)
    fgnSubNwVBox$packStart(fgnSubNwFrame, expand=FALSE,10)
    tabSubNetwork$packStart(fgnSubNwVBox, expand=FALSE,10)

    
    return(list(tabSubNetwork=tabSubNetwork))
}
