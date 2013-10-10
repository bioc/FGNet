getTerms <- function(results)
{	
	if(!is.list(results)) stop("Please provide the whole output from getResults_david() or getResults_gtLinker().")
	if("clusters" %in% names(results))
	{
		mgTerms <- getMGTerms(results$clusters, grType="Cluster")
	}
	if("metagroups" %in% names(results))
	{
		mgTerms <- getMGTerms(results$metagroups, grType="Metagroup")
	}
	
	terms <- lapply(mgTerms$termsDescriptions, function(x) 
		{
		temp <- rbind(x)[,"Description", drop=FALSE]
		colnames(temp) <- "Terms"
		rownames(temp) <- NULL
		return(temp)
		})
	return(terms)	
}
