\name{FGNet_GUI}
\alias{FGNet_GUI}
\title{
FGNet graphical user interface
}
\description{
Provides a graphical user interface (GUI) to most FGNet functionalities.
}
\usage{
FGNet_GUI(geneList = NULL)
}
\arguments{
  \item{geneList}{
vector. If provided, assigns the value to the genes field. It can be a character vector containing the gene list, or a named numeric vector with the gene expression.
}
}
\details{
To generate the functional network, first execute or import the Functional Analisis of the gene list with one of the tools in Tab "1 - FEA", then generate the network or the report in Tab "2 - Network".
}
\value{
Opens the GUI. No value is returned. The results of the analyses will be saved in the current working directory.
}

\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{
\dontrun{

FGNet_GUI()

# To directly input a gene list (i.e. from a previous analysis):
geneList <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W", 
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W",
    "YOL149W", "YOR249C")
# Optional gene expression
geneExpr <- setNames(c(rep(1,10),rep(-1,5)), geneList)

FGNet_GUI(geneExpr)
}
}