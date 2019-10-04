\name{FSARMA.fit}
\alias{FSARMA.fit}
\title{
Fit Seasonal ARMA
}

\description{
Fit basic seaonal ARMA model (more general than regular SARMA).
}

\usage{
FSARMA.fit(X,U,p,q,S1,S2,W=NULL,k=2)
}

\arguments{
\item{X}{data vector or a n by 1 matrix}

\item{U}{n-burnin, upper bound for seasonality}

\item{p}{AR lags}

\item{q}{MA lags}

\item{S1}{AR seasonality set}

\item{S2}{MA seasonality set}

\item{W}{exogenous matrix p by n}

\item{k}{IC penalty. k=2 is AIC; k=-1 select by CV}
}

\details{
The model is estimated by BFGS algorithm in optim(). Note that in univariate ARMA estimation, quasi-Newton method usually provide a robust result rather than aggressive ML with second order algorithms.

The algorithm optimize conditional likelihood based on burned in samples. This is specified by argument U. U has to be greater than p,q or any element in seasonality terms.

For models that the algorithm doesn't converge, the aic value is recorded as Inf.
}

\value{

\item{U}{n-burnin}

\item{p}{AR lags}

\item{phi}{estimated coefficients of AR lags terms}

\item{q}{MA lags}

\item{psi}{estimated coefficients of MA lags terms}

\item{r1}{length of S1}

\item{S1}{AR seasonality set}

\item{tau1}{estimated coefficients of AR seasonality set}

\item{r2}{length of S2}

\item{S2}{MA seasonality set}

\item{tau2}{estimated coefficients of MA seasonality set}

\item{gamma}{estimated coefficients of exogenous variables}

\item{sigma}{estimated sigma of white noise}

\item{aic}{if k=2, AIC value; if k=-1, CV error}

}

\author{
Tianyang Xie
}
