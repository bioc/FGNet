
# Returns color by order of metagroup (first color = first metagroup in list...)
setColors <- function(namesMg, transparency=0)
{
    nMetagroups <- length(namesMg)
    mgAlpha <- 1-transparency
    if(nMetagroups <5)
    {    
        colores <- c("#FF0000","#FFFF00","#0000FF", "#FF00E6")
        mgAlpha <- as.hexmode(round(mgAlpha*255))
        colores <- paste(colores, mgAlpha, sep="")
        
        colores <- colores[1:nMetagroups]
    }else
    {        
        secondStart <- 0.55
        if(nMetagroups > 8) secondStart <- 0.5

        n1 <- trunc(nMetagroups/2)
        n2 <- length((trunc(nMetagroups/2)+1):nMetagroups)
        
        colores <- rainbow(n1, end=0.3, s=1, v=1, alpha=mgAlpha)
        colores <- c(colores, rainbow(n2, start=secondStart, end=0.85, s=1, v=1, alpha=mgAlpha))
    }
    names(colores) <- namesMg
    
    return(colores)
}

#barplot(rep(1,length(colores)), col=as.character(colores), axes=FALSE)
