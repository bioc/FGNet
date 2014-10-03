## Retrieves the results form GeneTermLinker analysis
# Input: jobID ...
# Output: Results data.frame list

test_fea_gtLinker_getResults<- function()
{
	# Check the vignette/examples jobID still exists in the main GtLinker server and the function works propperly
	result <- fea_gtLinker_getResults(jobID=1639610)
	checkTrue(is.list(result))
	checkEquals(length(result), 4)
	checkTrue(all(names(result) %in% c("queryArgs","metagroups", "geneTermSets", "fileName")))
	checkTrue(is.data.frame(result$metagroups))

	# Check wrong ID
	# checkTrue(is.null(getResults_gtLinker(jobID=193283,serverWeb="http://gtlinker.cnb.csic.es")))
	
	# Check wrong URL
	# checkException(getResults_gtLinker(jobID=193283,serverWeb="http://gtlinker.cnb.csic2.es"), "Either the web server URL is wrong or the server is down")
}
