
clDescriptions <- function(terms, org=NULL) 
{
	terms <- sapply(terms, function(x) strsplit(x, split=":"))
	
	# Get term description		
	descripciones <- t(sapply(terms, function(x)
	{ 
		# Annotation
		annot <- x[1]
		
		# Description (Term name)
		descr <- capitalize(x[length(x)])
		
		# GO
		goID <- NA
		if(x[1] == "GO") goID <- x[2]
		
		# Link (Only for Kegg and Interpro)
		link <- NA
		# Interpro
		# David: 			IPR000719:Protein kinase, core
		# gtLinker: 	IPR000980:SH2 motif
		if(substring(x[1], 1, 3) == "IPR")
		{
			link <- paste("http://www.ebi.ac.uk/interpro/entry/", x[1], sep="")
			annot <- "InterPro"
		}
					
		# Kegg
		# David: 			K:xtr04330:Notch signaling pathway 
		if(x[1] == "KEGG") link <- paste("http://www.genome.jp/kegg-bin/show_pathway?org_name=", substring(x[2], 1, 3), "&mapno=", substring(x[2], 4, nchar(x[2])), sep="") 
		
		# gtLinker: 	Kegg:03040:Spliceosome 
		if(x[1] == "Kegg") link <- paste("http://www.genome.jp/kegg-bin/show_pathway?org_name=", org, "&mapno=", x[2], sep="")
		
		# SMART (Appears in david even if it was not explicitly requested)
		# David:			SM00181:EGF
		if(substring(x[1], 1, 3) == "SM0") link <- paste("http://smart.embl.de/smart/do_annotation.pl?DOMAIN=", x[1], sep="")			
		
		### NUEVOS
		if(length(grep("REACT", x[1]))>0) link <- paste("http://www.reactome.org/cgi-bin/link?SOURCE=Reactome&ID=", sub("REACT_", "REACT:", x[1]), sep="")	
		
		
		return (c(annot, descr, link, goID))
	}))
	colnames(descripciones) <- c("Annotation", "Description", "Link", "goID")
	return(descripciones)
}