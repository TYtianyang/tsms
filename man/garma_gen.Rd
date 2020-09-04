\name{garma.gen}
\alias{garma.gen}

\title{
Generate from generalized ARMA
}

\description{
Generate a sample series from a specified generalized ARMA model.
}

\usage{
garma.gen(n, U, p, q, S1 = NULL, S2 = NULL, phi = NULL, psi = NULL, tau1 = NULL, tau2 = NULL, noise = 1)
}

\arguments{

\item{n}{length of series}

\item{U}{number of data points burned-in}

\item{p}{number of regular AR components}

\item{q}{number of regular MA components}

\item{S1}{a set of lag orders of additional AR components}

\item{S2}{a set of lag orders of additional MA components}

\item{phi}{coefficients of regular AR components}

\item{psi}{coefficients of regular MA components}

\item{tau1}{coefficients of additional AR components}

\item{tau2}{coefficients of additional MA components}

\item{noise}{standard deviation of white noise}
}

\details{
Please note that the generator will automatically filter the generated series: those diverge dramatically will be ignored. However, it still does not gaurantee stationarity.
}

\value{
\item{X}{generated series in a n by 1 matrix}
  
\item{p}{number of regular AR components}

\item{q}{number of regular MA components}

\item{S1}{a set of lag orders of additional AR components}

\item{S2}{a set of lag orders of additional MA components}
  
\item{beta}{coefficients vector}
}

\author{
Tianyang Xie
}
