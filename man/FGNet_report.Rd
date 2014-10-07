\name{FGNet_report}
\alias{FGNet_report}
\alias{report_gtLinker}
\alias{report_david}
\title{
FGNet report
}
\description{
Generates an HTML report with several views of the Functional Network and complementary analyses.
}
\usage{
FGNet_report(feaResults, geneExpr = NULL, plotExpression = "border", 
onlyGoLeaves = TRUE, plotGoTree = TRUE, plotKeggPw = TRUE, 
filterAttribute = NULL, filterOperator = "<", filterThreshold = 0)
}
\arguments{
  \item{feaResults}{
list or data.frame/matrix. Output from one of the FEA functions.
}
  \item{geneExpr}{
numeric. Named vector with the relative expression value of the gene (node). 0 is taken as reference, positive values will be plotted red, negative values green.
}
  \item{plotExpression}{
character. Determines the way to plot the expression: "border" adds a red or green border to the node, "fill" colors the whole with the expression color instead of the metagroup color.
}
  \item{onlyGoLeaves}{
logical. If TRUE only terminal GO terms (leaves in the ontology tree) will be included in the cluster list.
}
  \item{plotGoTree}{
logical. If TRUE plots containing the terms in their position within the GO ontology (tree) will be generated.
}
  \item{plotKeggPw}{
logical. If TRUE local plot for KEGG pathway terms will be generated. If FALSE KEGG pathway terms will link to KEGG website.
In the current version local plots are only available for ENTREZ gene IDs.
}
  \item{filterAttribute}{
character or data.frame. Attribute to filter the clusters/metagroups. Filtered clusters/metagroups will not be included in the matrices (and subsequent networks).
Its value should be the data.frame column to use for filtering. It can be provided as character (column name) or data.frame (subset of the data.frame with drop=FALSE).
}
  \item{filterOperator}{
character. Logical operator used for filtering. i.e. ">" (bigger than), "<=" (smaller or equal than), "==" (equal), "!=" (different), "\%\%" (included in),... 
The evaluation order is left to right: filterAttribute ">" filterThreshold, will filter out clusters with filter attribute bigger than the threshold.
}
  \item{filterThreshold}{
numeric. Sets the value to compare to.
}
}
\value{
Generates the HTML report in the current directory.
}
\section{Warning}{
Reactome ID change depending on the database version. Links to reactome website are created to ease the analysis, but in case the version used for the enrichment and the website's do not match, broken links or misleading links might appear.  
}

\seealso{
Functional enrichment analysis functions:
\itemize{
    \item{\code{\link[FGNet:fea_david]{fea_david()}} (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gtLinker]{fea_gtLinker()}} & \code{fea_gtLinker_getResults()}
          (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gage]{fea_gage()}}}
    \item{\code{\link[FGNet:fea_topGO]{fea_topGO()}} 
    
   To import results from a previous/external FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.}
}


Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{
\dontrun{

# Report with diferent tools:

##########################
# DAVID & TopGO
geneList <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W", 
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W", 
    "YOL149W", "YOR249C")
geneExpr <- setNames(c(rep(1,10),rep(-1,5)), geneList) 

library(org.Sc.sgd.db)
geneLabels <- unlist(as.list(org.Sc.sgdGENENAME)[geneList])
names(geneExpr) <- geneLabels[names(geneExpr)] 

# DAVID
results_David <- fea_david(geneList, geneLabels=geneLabels, email="example@email.com")
FGNet_report(results_David, geneExpr=geneExpr) 

# TopGO
results_topGO <- fea_topGO(geneList, geneIdType="ENSEMBL", 
    geneLabels=geneLabels, organism="Sc") 
FGNet_report(results_topGO, geneExpr=geneExpr)   

##########################
# Gage
library(gage); data(gse16873)
results_gage <- fea_gage(eset=gse16873, 
    refSamples=grep('HN',colnames(gse16873), ignore.case =T), 
    compSamples=grep('DCIS',colnames(gse16873), ignore.case=T), 
    geneIdType="ENTREZID", organism="Hs", annotations="REACTOME")
FGNet_report(results_gage)

##########################
# Gene-Term Linker: 

# Execute new query:
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
     "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
     "GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
     "LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
     "PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
     "RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
     "SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")
# Optional expression (1=UP, -1=DW):
genesYeastExpr <- setNames(c(rep(1,29), rep(-1,30)), genesYeast)


jobID <- fea_gtLinker(geneList=genesYeast,organism="Sc")

# Load existing query:
jobID <- 3907019

results_gtLinker <- fea_gtLinker_getResults(jobID=jobID)
FGNet_report(results_gtLinker, geneExpr=genesYeastExpr)


}
}