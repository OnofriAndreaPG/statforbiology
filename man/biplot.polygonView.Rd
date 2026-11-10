\name{biplot.polygonView}
\alias{biplot.polygonView}

\title{
Draws a biplot after GGE analysis and superimposes a 'which-won-where' pattern, as defined in Yan and Tinker (2006).
}
\description{
This function is reserved for the class 'GGEobject' as obtained from the 'GGE()' function in the 'statforbiology' package.
}
\usage{
biplot.polygonView(obj, xlim=NULL, ylim=NULL, elabels=NULL,
               glabels=NULL, quad=FALSE, cexG=0.9,
               cexE=0.9,
               xlab=NULL, ylab=NULL, font=1, ...)
}
\arguments{
  \item{obj}{a 'GGEobject' objects}
  \item{xlim}{graphical parameter as in the 'plot' method}
  \item{ylim}{graphical parameter as in the 'plot' method}
  \item{elabels}{labels for the environments}
  \item{glabels}{labels for the genotypes}
  \item{quad}{logical. If TRUE, plots the axes}
  \item{cexG}{graphical parameter: the 'cex' for the genotype labels}
  \item{cexE}{graphical parameter: the 'cex' for the environment labels}
  \item{xlab}{graphical parameter as in the 'plot' method}
  \item{ylab}{graphical parameter as in the 'plot' method}
  \item{font}{graphical parameter as in the 'plot' method. Relating to the genotype labels}
  \item{...}{Other additional arguments}
}
\value{No value is returned}

\author{
Andrea Onofri
}
\examples{
dataset <- getAgroData("WinterWheat")
tab <- with(dataset, GGE(Yield, Genotype, Year, Block, PC = 2))
biplot.polygonView(tab, xlim = c(-2, 2), ylim = c(-2,2))
}
