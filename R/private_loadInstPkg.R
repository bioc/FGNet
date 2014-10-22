
# usage:
# if(!loadInstPkg("gage")) stop("Package ... is not available.")
# if(!loadInstPkg("gage", parentWindow=argsList$parentWindow)) stop("Package gage is not available.")
# 
loadInstPkg <- function(pkgName, parentWindow=NULL)
{
    ret <- FALSE
    
    if(!pkgName %in% rownames(installed.packages()))
    {
        #if(interactive() && any(grepl("RGtk2", unlist(sessionInfo()))))
        if(!is.null(parentWindow))
        {
            # GUI dialog
            dialog <- gtkDialogNewWithButtons("Package installation required", parentWindow,
                                              c("modal", "destroy-with-parent"), 
                                              "gtk-yes", GtkResponseType["yes"], 
                                              "gtk-no", GtkResponseType["no"],
                                              show=TRUE)
            dialog[["vbox"]]$add(gtkLabel(paste("Install package '", pkgName,"'?", sep="")))
            response <- dialog$run() 
            dialog$destroy()
            
            if(response == GtkResponseType["yes"]) response <- "y"
        }else
        { 
            # Console
            response <- readline(paste("Package '",pkgName, "' is required but not installed. Install now? (y/n) \n", sep=""))
        }
        
        if(grepl("y", tolower(response)))
        {
            tryCatch({
                source("http://bioconductor.org/biocLite.R")
                biocLite(pkgName)
            }, error = function(e) {
                stop(paste("It is not possible to install the package. Check your internet connection.", e, sep="\n"))
            })
            ret <- TRUE
        }
    }
    
    # pos = "package:base"  # to "mask as little as possible"
    ret <- suppressWarnings(library(pkgName, character.only=TRUE, logical.return=TRUE, quietly=TRUE, pos = "package:base")) # TRUE/FALSE
    
    return(ret)
}
