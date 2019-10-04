\name{HSARMA.auto}
\alias{HSARMA.auto}

\title{
Hybrid Seasonal ARMA
}

\description{
Detect potential seasonality periods by spectrum analysis, fit models on blurred seasonality periods terms, select best model, and predict. For centered and detrended univariate time series.
}

\usage{
HSARMA.auto(X,W=NULL,U,p_range,q_range,r=1,S=NULL,
                       blur.out=c(2,2),sar=T,sma=F,sfourier=T,order=1,pred_t=0,
                       k=2,level=0.05,multicore=T)
}

\arguments{
\item{X}{data vector or a n by 1 matrix}

\item{W}{exogenous matrix p by (n+pred_t)}

\item{U}{n-burnin, upper bound for seasonality}

\item{p_range}{range of p, AR lags}

\item{q_range}{range of q, MA lags}

\item{r}{# seasonality periods}

\item{S}{Pre-specified seasonal terms;if NULL, suggested by frequency domain techniques}

\item{blur.out}{2-elements vector, left and right bounds of seasonality intervals}

\item{sar}{logical, including seasonal AR components or not}

\item{sma}{logical, including seasonal MA components or not}

\item{sfourier}{logical, including triangular components or not}

\item{order}{order of triangular components}

\item{pred_t}{n-predict out of samples. =NULL then will not predict}

\item{k}{selection criteria. k=2 is AIC; k=-1 is CV}

\item{level}{level of z-score test on spectral spikes}

\item{multicore}{logical, whether to employ cluster}

}

\details{
The model is estimated by BFGS algorithm in optim(). Note that in univariate ARMA estimation, quasi-Newton method usually provide a robust result rather than aggressive ML with second order algorithms.

The algorithm optimize conditional likelihood based on burned in samples. This is specified by argument U. U has to be greater than p,q or any element in seasonality terms.

For models that the algorithm doesn't converge, the aic value is recorded as Inf.

The function make static prediction within each observation period, while use dynamic prediction for out-sampled period.
}

\value{

\item{obj}{best fitted model object}

\item{aic_panel}{selection panel for all parameters' combination}

\item{prediction}{a vector recording predicted values. Length= n + pred_t}

}

\author{
Tianyang Xie
}


