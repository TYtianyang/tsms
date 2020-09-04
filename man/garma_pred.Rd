\name{garma.pred}
\alias{garma.pred}

\title{
Predict generalized ARMA
}

\description{
Predict the fitted generalized ARMA model.
}

\usage{
garma.pred(obj,X,W=NULL,pred_t)
}

\arguments{

\item{obj}{the fitted ARMA model object}

\item{X}{a data vector or n by 1 matrix}

\item{W}{exogenous variable matrix p by (n+pred_t) matrix}

\item{pred_t}{n-outsample prediction length}
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



