# Used in GUI and REPORT
# if(is.null(jobName)) jobName <- getJobName(results$fileName) 

getJobName <- function(fileName)
{    
    if(is.numeric(fileName)) # jobID?
    {
        jobName <- as.character(fileName)
    }else
    {
        # Remove folder/location
        jobName <- strsplit(fileName, .Platform$file.sep, fixed=TRUE)[[1]]
        jobName <- jobName[length(jobName)]
        
        # Gtlinker
        jobName <- gsub("_formatted_metagroups.txt","", jobName)
        jobName <- gsub("_formatted_gtSets.txt","", jobName)
        jobName <- gsub("_global_overview.txt","", jobName)
        
        # Other
        jobName <- gsub("_formatted.txt","", jobName)
        jobName <- gsub(".txt","", jobName)
        jobName <- gsub("_feaResults.RData","", jobName)
        jobName <- gsub(".RData","", jobName)
        
    }
    return(jobName)
}
