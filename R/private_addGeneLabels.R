

addGeneLabels <- function(tablaGeneTermSets, geneLabels)
{
    if(!is.null(geneLabels))
    {
        # Shouldnt happen...
        if(length(geneLabels) != length(unique(names(geneLabels)))) stop("geneLabels IDs are not unique.")
     
        # Reorder colums
        colGenes <- which(colnames(tablaGeneTermSets) == "Genes")
        colnames(tablaGeneTermSets)[colGenes] <- "GenesIDs"
        tablaGeneTermSets <- tablaGeneTermSets[,c(colnames(tablaGeneTermSets)[-colGenes], "GenesIDs")] 
        
        # Replace gene IDs in temp var
        tmpGenes <- sapply(as.character(tablaGeneTermSets[,"GenesIDs"]), function(clGenes){
            paste(geneLabels[strsplit(as.character(clGenes),split=",")[[1]]], collapse=",")    
        }) 
        
        # add
        tablaGeneTermSets <- cbind(tablaGeneTermSets, Genes=tmpGenes)
    }
    return(tablaGeneTermSets)
}
