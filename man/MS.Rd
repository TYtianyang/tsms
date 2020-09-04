\name{MS}
\alias{MS}

\title{
Multiple Seasonality modeling procedure
}

\description{
Detect potential seasonality periods by spectrum analysis, fit models on blurred seasonality periods terms, select best model, and predict. Now only for centered and detrended univariate time series.
}

\usage{
MS(X, W = NULL, U, p_range, q_range, r = 1, S = NULL, 
    blur.out = c(2, 2), sar = T, sma = T, sfourier = F, order = 1, 
    sfactor = F, pred_t = 0, crit = "BC", level = 0.05, 
    multicore = T)
}

\arguments{
\item{X}{a data vector or a n by 1 matrix}

\item{W}{exogenous matrix p by (n+pred_t)}

\item{U}{number of data points burned-in. Upper bound for seasonality, i.e. U > max(S1,S2)}

\item{p_range}{a range of regular AR lag to be selected}

\item{q_range}{a range of regular MA lag to be selected}

\item{r}{number of seasonality periods allowed}

\item{S}{pre-specified seasonal ARMA components;if NULL, suggested by information criterion and MS procedure will be performed}

\item{blur.out}{2-elements vector, left and right bounds of blurring seasonality intervals}

\item{sar}{logical, including seasonal AR components or not}

\item{sma}{logical, including seasonal MA components or not}

\item{sfourier}{logical, including triangular components or not}

\item{order}{order of triangular components}
  
\item{sfactor}{logical, including binary seasonal components or not}

\item{pred_t}{n-predict out of samples. if NULL then will not predict}

\item{crit}{selection criterion. crit = c('BC','AIC','BIC')}

\item{level}{level of z-score test on spectral spikes. Used to control significance of potential seasonality in the first step}

\item{multicore}{logical, whether to employ parallel computation}

}

\details{
The model is estimated by BFGS algorithm in optim(). Note that in univariate ARMA estimation, quasi-Newton method usually provide a robust result rather than aggressive ML with second order algorithms.

The algorithm optimize conditional likelihood based on burned in samples. This is specified by argument U. U has to be greater than p,q or any element in seasonality terms.

For models that have diverging estimation, the aic value will be recorded as Inf.

The function make static prediction within each observation period, while use dynamic prediction for out-sampled period.
}

\value{

\item{obj}{best fitted model object}

\item{ic_panel}{selection panel for all parameters combination}

\item{prediction}{a vector recording predicted values. Length= n + pred_t}

}

\author{
Tianyang Xie
}


