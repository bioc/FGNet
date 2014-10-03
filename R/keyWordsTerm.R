# Used by createHTML and clustersTable
# nChar number of characters to show (trims the result)
# i.e. termsDescriptions <- getTerms(results)
keywordsTerm <- function(termsDescriptions, nChar=30)
{    
    if(!is.list(termsDescriptions)) termsDescriptions <-list(termsDescriptions)
    if(ncol(termsDescriptions[[1]]) > 1) termsDescriptions <- lapply(termsDescriptions, function(x) cbind(x[,"TermDescription", drop=FALSE]))
    
    # Search for a representative term (based only on key words! not function!)
    reprTerms <- rep(NA, length(termsDescriptions))
    for(index in 1:length(termsDescriptions))
    {        
        termsMg <- termsDescriptions[[index]]
        
        rownames(termsMg) <- NULL
        splittedTerms <- unlist(strsplit(termsMg, " "))
        #splittedTerms <- unlist(strsplit(splittedTerms, "-", fixed=TRUE))
        splittedTerms <- unlist(strsplit(splittedTerms, "/", fixed=TRUE))
        splittedTerms <- unlist(strsplit(splittedTerms, ",", fixed=TRUE))
        splittedTerms <- unlist(strsplit(splittedTerms, ".", fixed=TRUE))
        splittedTerms <- unlist(strsplit(splittedTerms, "(", fixed=TRUE))
        splittedTerms <- unlist(strsplit(splittedTerms, ")", fixed=TRUE))
        
        sortedFreq <- sort(table(tolower(splittedTerms)), decreasing=TRUE)
        
        commonWords <- sortedFreq[sortedFreq >= quantile(sortedFreq, 0.9)]
        commonWords <- commonWords[grepl("[[:alnum:]]", names(commonWords))]
        
        selectedTerm <- NULL
        if(any(tolower(termsMg) %in% names(commonWords)))
        {
            selectedTerm <- names(commonWords)[which(names(commonWords) %in% tolower(termsMg))[1]] # Select the one with highest freq.
            
            if(commonWords[1] >= commonWords[selectedTerm]*2) selectedTerm <- NULL     
            selectedTerm <- capitalize(selectedTerm)
        }
        if(is.null(selectedTerm))
        {
            termsWithCommonWords <- sapply(names(commonWords), function(x) grep(x, tolower(termsMg)))
            countTermsWords <- table(unlist(termsWithCommonWords)) #Number of common words per term        
            selectedTerms <- termsMg[as.numeric(names(countTermsWords[countTermsWords == max(countTermsWords)])),] # Terms with most common words
            ncharTerms <- nchar(selectedTerms)
            selectedTerm <- selectedTerms[ncharTerms==min(ncharTerms)][1]
        }
        
        reprTerms[index] <- selectedTerm
        if(nchar(reprTerms[index])>nChar) reprTerms[index] <- paste(substr(reprTerms[index], 1, nChar), "-", sep="")
    }
    if(!is.null(names(termsDescriptions))) names(reprTerms) <- names(termsDescriptions)
    return(reprTerms)
}