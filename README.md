# The 'tsms' R package
The tsms R package introduces a tool box for generalized ARMA model, as well as the MS modeling procedure for univariate time series with multiple seasonality. 

## Getting Started

First install the devtools package

install.packages("devtools")

library("devtools")

Then install this package

install_github('TYtianyang/tsms')

Note: For R version > 4.0.0 user, avoid auto-converting warnings to error issues from older packages (e.g. gridExtra) by the following command before installation

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

## Using This Package

To see the available function to use, type 

ls("package:tsms")
