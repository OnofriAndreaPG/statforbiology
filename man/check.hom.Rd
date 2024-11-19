\name{check.hom}
\alias{check.hom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check linear models for homoscedasticity
}
\description{
This function takes a linear model object as an argument and checks whether the residuals are homoscedastic, in relation to a stratification variable or covariate, that is given as an argument.
}
\usage{
check.hom(obj, var, alternative)
}
\arguments{
  \item{obj}{a linear model object fitted with lm()}
  \item{var}{a vector containing a stratification variable. If missing, the fitted values from obj are taken as a covariate}
  \item{alternative}{The null (homogeneous variances) is tested against the alternative, that is (i) different variances for each level of the stratification variable given as argument ("varIdent"), (ii) variance is a power function of a covariate, that is given as argument ("varPower"; the fitted values are taken as the covariate, in case the argument 'var' is missing and (iii) the variance is an exponential function of a covariate, that is given as argument ("varExp"; the fitted values are taken as the covariate, in case the argument 'var' is missing)}
}
\details{The function fits a gls with same structure as the input model, together with a heteroscedastic gls, where residuals are allowed to be heteroscedastic, according to the 'alternative' argument. The two models are compared by a LRT}
\value{
\item{LRT}{the value of LRT}
\item{LRT}{the P-value of LRT}
\item{aovtable}{a summary table for the LRT}
\item{modHet}{the gls object containing the heteroscedastic fit}
}
\author{
Andrea Onofri
}

\examples{
fileName <- "https://www.casaonofri.it/_datasets/FGP_rape.csv"
library(statforbiology)
dataset <- read.csv(fileName)
dataset[,1:5] <- lapply(dataset[,1:5], factor)
mod <- lm(FGP ~ Genotype * Run, data = dataset)
check <- check.hom(mod, Run)
check$aovtable
}
