\name{FSARMA.pred}
\alias{FSARMA.pred}

\title{
Predict on fitted Seasonal ARMA
}

\description{
Predict on fitted Seasonal ARMA
}

\usage{
FSARMA.pred(obj,X,W=NULL,pred_t)
}

\arguments{

\item{obj}{best fitted FSARMA.fit object}

\item{X}{data vector or n by 1 matrix}

\item{W}{exogenous matrix p by (n+pred_t) matrix}

\item{pred_t}{n-outsample prediction}
}

\details{
The function make static prediction within observation period, while use dynamic prediction for out sampled period. For burned-in observations, the predicted value will be recorded as 0.
}

\value{

\item{pred_X}{a vector recording predicted values. Length = n + pred_t}
}

\author{
Tianyang Xie
}



