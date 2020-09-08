# The 'tsms' R package
The tsms R package introduces a tool box for generalized ARMA model, as well as the MS modeling procedure for univariate time series with multiple seasonality. The MS modeling procedure is a simple while well performed modeling machine for time series with multiple seasonality. It detects potential seasonality orders by Discrete Fourier Transformation (DFT), fits multiple generalized ARMA model parallelly by the BFGS algorithm, and selects the best model by minimizing the information criterion (IC). 

Unlike most of the methods for multiple seasonality, MS doesn't require pre-specification of the seasonality orders (e.g. Weekly, monthly or yearly seasonlity cycles). It's 'automatic' in a way that the user only need to specify the number of seasonlity cycles. From simulation and empirical studies, we claim that MS has superior performance than the benchmark methods (e.g. ARIMA, Facebook Prophet, TBATS), and comparable performance to SOTA methods (e.g. LSTM). 

More details can be referred to https://arxiv.org/abs/2008.12340. Feel free to contact the authors: Tianyang Xie at xie00039@umn.edu, Jie Ding at dingj@umn.edu .

## Getting Started

install.packages("devtools")

library("devtools")

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

install_github('TYtianyang/tsms')

Note: For R version > 4.0.0 user, we avoid auto-converting warnings to error issues from older packages (e.g. gridExtra) by the the command: 

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

## Using This Package

To see the available function to use, type 

ls("package:tsms")
