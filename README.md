# The 'tsms' R package
The tsms R package introduces a tool box for generalized ARMA model, as well as the MS modeling procedure for univariate time series with multiple seasonality. The MS modeling procedure is a simple while well performed modeling machine for time series with multiple seasonality. It detects potential seasonality orders by Discrete Fourier Transformation (DFT), fits multiple generalized ARMA model parallelly by the BFGS algorithm, and selects the best model by minimizing the information criterion (IC). 

Unlike most of the methods for multiple seasonality, MS doesn't require pre-specification of the seasonality orders (e.g. Weekly, monthly or yearly seasonlity cycles). It's 'automatic' in a way that the user only need to specify the number of seasonlity cycles. From simulation and empirical studies, we claim that MS has superior performance than the benchmark methods (e.g. ARIMA, Facebook Prophet, TBATS), and comparable performance to SOTA methods (e.g. LSTM). 

More details can be referred to https://arxiv.org/abs/2008.12340. Feel free to contact the authors: Tianyang Xie at xie00039@umn.edu, Jie Ding at dingj@umn.edu .

## Getting Started
```{r}
install.packages("devtools")

library("devtools")

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

install_github('TYtianyang/tsms')
```


__Note: For R version > 4.0.0 user, we avoid auto-converting warnings to error issues from older packages (e.g. gridExtra) by the the command:__ 

__Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")__

## Using This Package

To see the available function to use, type 

```{r}
ls("package:tsms")
```

For explicit usage of each function, please refer to [Explore Generalized ARMA with TSMS](https://github.com/TYtianyang/tsms/blob/master/vignettes/explore_generalized_arma_with_tsms.pdf)

For a MS procedure real-data example, please check [A MS Example](https://github.com/TYtianyang/tsms/blob/master/vignettes/ms_example.pdf)
