\name{garma.fit}
\alias{garma.fit}
\title{
Fit generalized ARMA
}

\description{
Fit a generalized ARMA model for univariate time series. The ARMA is generalized in a way that AR and MA components can be specified in any lag orders.
}

\usage{
garma.fit(X,U,p,q,S1,S2,W=NULL,crit='AIC')
}

\arguments{
\item{X}{a data vector or a n by 1 matrix}

\item{U}{number of data points burned-in. Upper bound for seasonality, i.e. U > max(S1,S2)}

\item{p}{number of regular AR components}

\item{q}{number of regular MA components}

\item{S1}{a set of lag orders of additional AR components}

\item{S2}{a set of lag orders of additional MA components}

\item{W}{exogenous variable matrix p by n}

\item{crit}{selection criterion. crit = c('BC','AIC','BIC')}
}

\details{
The model is estimated by BFGS algorithm in optim(). Note that in univariate ARMA estimation, quasi-Newton method usually provide a robust result rather than aggressive ML with second order algorithms.

The algorithm optimize conditional likelihood based on burned in samples. This is specified by argument U. U has to be greater than p,q or any element in seasonality terms.

For models that have diverging estimation, the aic value will be recorded as Inf.
}

\value{

\item{U}{n-burnin}

\item{p}{number of regular AR components}

\item{phi}{estimated coefficients of regular AR components}

\item{q}{number of regular MA components}

\item{psi}{estimated coefficients of regular MA components}

\item{r1}{length of S1}

\item{S1}{a set of lag orders of additional AR components}

\item{tau1}{estimated coefficients of additional AR components}

\item{r2}{length of S2}

\item{S2}{a set of lag orders of additional MA components}

\item{tau2}{estimated coefficients of additional MA components}

\item{gamma}{estimated coefficients of exogenous variables}

\item{sigma}{estimated sigma of white noise}

\item{ic}{information criterion}

}

\author{
Tianyang Xie
}
