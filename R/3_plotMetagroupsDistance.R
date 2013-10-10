
# Distances: Binary
# Clustering: Single
plotMetagroupsDistance <- function(metagroupGenesMatrix)
{	
	mxDistancias <- NULL
	if(dim(metagroupGenesMatrix)[2]>2)	
	{
		distancias <- dist(t(metagroupGenesMatrix), method="binary")
		mxDistancias <- as.matrix(distancias)	
		clustering <- hclust(distancias, method="single")
		
		# Prepare layout
		nf <- layout(cbind(c(1,2),c(0,3)), widths=c(9,1), heights=c(3,8), TRUE)
		#layout.show(nf)
		
		# Plot dendogram
		par(mar = c(0,3,1,1))
		plot(clustering,  ann=FALSE)
		title("Metagroups distance")
		
		
		#Recolocar:
		mxDistancias <- mxDistancias[clustering$order, rev(clustering$order)]

		# Plot "image"		
		minDist <- unique(sort(mxDistancias))[2] #[1] = 0		
		if(minDist > 0)
		{
			cero <- which(mxDistancias == 0)
			uno  <- which(mxDistancias == 1)
			
			mxPlot <- (mxDistancias-(minDist))+0.12			
			mxPlot <- mxPlot/(max(mxPlot)/0.9)
			mxPlot <- trunc(mxPlot*10)
			
			mxPlot[cero] <- 0
			mxPlot[uno] <- 10
			
		} else mxPlot <- mxDistancias

		par(mar = c(3,3,1,1))
		
		if (library(RColorBrewer, logical.return=TRUE))
		{
			colors <- rev(brewer.pal(11, "RdBu"))
		}else
		{
			colors <- heat.colors(11)
		}
		image(mxPlot, col=colors, axes=FALSE)
		box()
		axis(1, at=seq(0,1, by=1/(length(colnames(mxPlot))-1)), labels = rownames(mxPlot))
		axis(2, at=seq(0,1, by=1/(length(colnames(mxPlot))-1)), labels = colnames(mxPlot))
		
 		# Leyenda
 		levs <- unique(sort(mxPlot))
 		labs <- sapply(lapply(levs, function(x) which(mxPlot %in% x)), function(x) round(mean(mxDistancias[x]), digits=2))
 		par(mar = c(3,1,15,0.5))
 		image(t(as.matrix(levs)),  col=colors, axes=FALSE)
		box()
	  axis(2, at=seq(0,1, by=1/(length(levs)-1)), labels= labs, cex.axis=0.8)

 		#par(mar= c(5, 4, 4, 2) + 0.1)  #- reset to default
		mxDistancias <- mxDistancias[,rev(colnames(mxDistancias))]
	}
	return(mxDistancias)
}
