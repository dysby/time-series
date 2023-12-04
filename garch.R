library(latex2exp)
library(patchwork)
library(ggpubr)
library(fGarch)
library(rugarch)
library(fable)
library(feasts)
library(kableExtra)
require(FinTS)
library(tidyverse, quietly = TRUE)

# Stylized facts:
#   The sample mean of the data is close to zero whereas the sample
# variance is of the order 10−4 or smaller;

# Exceedances of high/low thresholds tend to occur in clusters
# (dependence in the tails);

# Heavy-tailed marginal distributions;

# The sample ACF use to be negligible at all lags (with a possible exception of the first lag).
# On the contrary, the sample ACF of the absolute values or the squares are
# different from zero for a large number of lags and stay almost constant and
# positive for large lags;

# Asymmetric responses in the volatility: volatility tends to rise
# higher in response to negative shocks as opposed to positive shocks.

# http://eclr.humanities.manchester.ac.uk/index.php/R_GARCH
# https://github.com/msperlin/GARCH-RAC/blob/master/fcts/garch_fcts.R

# In the variance equation, we must pay attention to the value of the and . We
# expect them to be positive and their sum should not be more than one. In our case, their
# sum equals to 0,98. If this is not the case and the sum of ARCH and GARCH parameter
# is higher than one, the model has explosive behavior -- will keep increasingly its value
# infinitely. If this occurs, it is very likely that there was an error in the estimation of the
# model, probably in the input series.


adf <- function(x, ...) {
  out <- tseries::adf.test(x, k = 1)
  c(adf_stat = unname(out$statistic), adf_pvalue = out$p.value)
}

do_arch_test <- function(x, max_lag = 5) {
  do_single_arch <- function(x, used_lag)  {
    test_out <- FinTS::ArchTest(x, lags = used_lag)

    res_out <- tibble(Lag = used_lag,
                      `LMStatistic` = test_out$statistic,
                      `pvalue` = test_out$p.value)
  }

  tab_out <- bind_rows(map(1:max_lag, .f = do_single_arch, x = x))

  return(tab_out)
}

# ---- LOAD DATA ----
load_logreturns_xls <- function(path) {
  readxl::read_xls(path = path, skip = 3) %>%
    select(Close) %>%
    ts(., frequency = 1) %>%
    log() %>%
    diff()
}

load_ts_xls <- function(path) {
  readxl::read_xls(path = path, skip = 3) %>%
    select(Close) %>%
    ts(., frequency = 1)
}


EDP <- load_ts_xls("data/EDP RENOVAVEISprice.xls")
GALP <- load_ts_xls("data/GALP ENERGIA-NOMprice.xls")
NOS <- load_ts_xls("data/NOSSGPSprice.xls")
NOVABASE <- load_ts_xls("data/NOVABASESGPSprice.xls")
MOTAENGIL <- load_ts_xls("data/MOTA ENGILprice.xls")


log_returnEDP <- load_logreturns_xls("data/EDP RENOVAVEISprice.xls"); colnames(log_returnEDP) <- c('log_returnEDP')
log_returnGALP <- load_logreturns_xls("data/GALP ENERGIA-NOMprice.xls"); colnames(log_returnGALP) <- c('log_returnGALP')
log_returnNOS <- load_logreturns_xls("data/NOSSGPSprice.xls"); colnames(log_returnNOS) <- c('log_returnNOS')
log_returnNOVABASE <- load_logreturns_xls("data/NOVABASESGPSprice.xls"); colnames(log_returnNOVABASE) <- c('log_returnNOVABASE')
log_returnMOTAENGIL <- load_logreturns_xls("data/MOTA ENGILprice.xls"); colnames(log_returnMOTAENGIL) <- c('log_returnMOTAENGIL')

# table stats

table.Stats(log_returnEDP)

ggplot(data=EDP, aes(x=x, y=y)) +
  geom_line()

log_returnEDP %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnGALP %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnNOS %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnNOVABASE %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnMOTAENGIL %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
# MOTAENGIL -Heteroscedasticity and bigger negative shocks (skewed)
# MOTAENGIL seems to have volatility clustering.
log_returnEDP^2 %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnGALP^2 %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnNOS^2 %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnNOVABASE^2 %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()
log_returnMOTAENGIL^2 %>% ggplot(aes(x=x, y=y)) + ylim(c(-0.10, 0.10)) + geom_line()

# Testing for ARCH/GARCH Effects
do_arch_test(log_returnEDP)
do_arch_test(log_returnGALP)
do_arch_test(log_returnNOS)
do_arch_test(log_returnNOVABASE)
do_arch_test(log_returnMOTAENGIL)

ArchTest(log_returnEDP)

Box.test(coredata(log_returnEDP^2), type="Ljung-Box", lag = 12)
Box.test(coredata(log_returnGALP^2), type="Ljung-Box", lag = 12)
Box.test(coredata(log_returnNOS^2), type="Ljung-Box", lag = 12)
Box.test(coredata(log_returnNOVABASE^2), type="Ljung-Box", lag = 12)
Box.test(coredata(log_returnMOTAENGIL^2), type="Ljung-Box", lag = 12)
# MOTAENGIL have ARCH effect

as_tsibble(NOS) %>% features(value, list(unitroot_kpss, adf))
as_tsibble(log_returnNOS) %>% features(value, list(unitroot_kpss, adf))

# series_tb <- tibble(EDP = log_returnEDP,
#                     GALP = log_returnGALP,
#                     NOS = log_returnNOS,
#                     #NOVABASE = log_returnNOVABASE,
#                     MOTAENGIL = log_returnMOTAENGIL)

log_returns_list <- list(log_returnEDP, log_returnGALP, log_returnNOS, log_returnNOVABASE, log_returnMOTAENGIL)
log_returns_names <- c("EDP", "GALP", "NOS", "NOVABASE", "MOTAENGIL")

series_stats <- lapply(log_returns_list,
                       function(lg) {
                         FinTS::FinTS.stats(lg)
                       }
) %>%
  bind_rows() %>%
  cbind(log_returns_names) %>%
  rename("Series (log_return)" = log_returns_names) %>%
  select(9, 2:8) %>%
  mutate_at(vars("Mean", "Standard.Deviation", "Skewness", "Excess.Kurtosis", "Minimum", "Maximum"),
            list(~ round(., digits = 4)))

# series_tb %>%
#   summarise_all(function(lg) {
#     FinTS::FinTS.stats(lg)
# })

log_returnEDP %>% forecast::Acf(main="Log returns", )
log_returnEDP %>% forecast::Acf(main="Log returns", type = "partial", ylim=c(-0.10, 1))

plots_acf_pacf <- function(s, main="Series") {
  par(mfrow = c(2, 2))
  s %>% forecast::Acf(main=main, ylim=c(-0.10, 1))
  s %>% forecast::Acf(main="", type = "partial", ylim=c(-0.10, 1))
}
plots_acf_pacf(log_returnEDP)
plots_acf_pacf(log_returnEDP^2)
# For lrEDP ACF and lrEDP^2 do not have signiticant lags so the series is uncorrelated,
# if a ARCH family model is to be fitted alpha1 coefficient will be nonsignificant because there is no dependence on the past values of squared series.
plots_acf_pacf(log_returnGALP)
plots_acf_pacf(log_returnGALP^2)
# Also for GALP
plots_acf_pacf(log_returnNOS)
plots_acf_pacf(log_returnNOS^2)
# Also for NOS
plots_acf_pacf(log_returnNOVABASE)
plots_acf_pacf(log_returnNOVABASE^2)
# Also for NOVABASE
plots_acf_pacf(log_returnMOTAENGIL)
plots_acf_pacf(log_returnMOTAENGIL^2)
# lrMOTAENGIL have some correlation with past squared values.

# !IMPORTANT
# see plot of ACF(log_returns) non-correlation, plot of ACF(log_returns^2) show correlations, therefore log_returns is not autocorrelated but is autodependent!
#
# https://cpb-us-w2.wpmucdn.com/blog.nus.edu.sg/dist/0/6796/files/2017/03/analysis-of-financial-time-series-copy-2ffgm3v.pdf
# Analysis of Financial Time Series, p100
# These two plots clearly suggest that the monthly returns are not serially
# independent. Combining the three plots, it seems that the returns are indeed serially
# uncorrelated, but dependent. Volatility models attempt to capture such dependence
# in the return series.

cbind(PerformanceAnalytics::table.Stats(log_returnEDP),
      PerformanceAnalytics::table.Stats(log_returnGALP),
      PerformanceAnalytics::table.Stats(log_returnNOS),
      PerformanceAnalytics::table.Stats(log_returnNOVABASE),
      PerformanceAnalytics::table.Stats(log_returnMOTAENGIL))

#                 log_returnEDP log_returnGALP log_returnNOS log_returnNOVABASE log_returnMOTAENGIL
# Observations         508.0000       508.0000      508.0000           507.0000            508.0000
# Arithmetic Mean       -0.0017        -0.0006       -0.0003            -0.0009             -0.0004
# Variance               0.0006         0.0006        0.0003             0.0003              0.0008
# Skewness              -0.0880        -1.0032       -0.6261             0.4085             -1.7197
# Kurtosis               1.5650         5.1624        8.0971             6.9793             23.2531

# From ACF analisys we identify ARMA(0,1) model, lets fit it to get the residuals and apply garch to the residuals.
m <- tseries::arma(log_returnEDP, order = c(7, 0), include.intercept = TRUE)
summary(m)

m <- tseries::arma(log_returnGALP, order = c(0, 1), include.intercept = FALSE)
summary(m)

# ---- ARCH ----
arch_fit <- garchFit(formula = . ~ garch(1,0), data = log_returnEDP, trace = FALSE)
show(arch_fit)
summary(arch_fit)

# garchFit(formula = . ~ garch(1, 0), data = log_returnEDP, trace = FALSE)
#   Estimate  Std. Error  t value Pr(>|t|)
# mu     -2.033e-03   1.036e-03   -1.961   0.0498 *
# omega   4.975e-04   4.494e-05   11.070   <2e-16 ***
# alpha1  1.167e-01   7.588e-02    1.537   0.1242
# Standardised Residuals Tests:
#   Statistic p-Value
# Jarque-Bera Test   R    Chi^2  53.15615  2.866041e-12
# Shapiro-Wilk Test  R    W      0.9854043 5.617446e-05
# Ljung-Box Test     R    Q(10)  16.67124  0.08196247
# Ljung-Box Test     R    Q(15)  20.15111  0.1662121
# Ljung-Box Test     R    Q(20)  23.79199  0.2515967
# Ljung-Box Test     R^2  Q(10)  19.12878  0.03865719
# Ljung-Box Test     R^2  Q(15)  28.22386  0.02020902
# Ljung-Box Test     R^2  Q(20)  30.121    0.06791668
# LM Arch Test       R    TR^2   6.472527  0.8904166
#
# Information Criterion Statistics:
#   AIC       BIC       SIC      HQIC
# -4.649665 -4.624682 -4.649734 -4.639868


# AUTO ARIMA -----


library(forecast)
armL=sapply(log_returns_list, function(x) unlist(auto.arima(x)$arma) )
pk=c(1,3,2) #<<< pick the rows to get AR/MA orders, row 1 is AR; 2 is MA
t(armL[pk,]) ### gives the auto.arima returned order

#       [,1] [,2] [,3]
# [1,]    0    0    0
# [2,]    0    0    0
# [3,]    0    0    0
# [4,]    0    0    0
# [5,]    1    0    1



# arch_fit@fit$ics # AIC and BIC
#plots_acf_pacf(arch_EDP@residuals, "ARCH(1) Residuals")

s <- log_returnEDP

arma.order <- c(1,1)
arma.order <- c(0,0)

arch.order <-  1:2
arch.names <-  paste("arch", arch.order, sep="")
fit.list <-  list()
for (p in arch.order) {
  arch.spec <-  ugarchspec(variance.model = list(garchOrder=c(p,0)),
                          mean.model = list(armaOrder=arma.order))
  fit.list[[p]] <- ugarchfit(spec=arch.spec, data=s)
}
names(fit.list) <-  arch.names

garch11.spec <-  ugarchspec(variance.model=list(garchOrder=c(1,1)),
                            mean.model = list(armaOrder=arma.order))
fit.list$garch11 <- ugarchfit(spec=garch11.spec, data=s)

garch11_std.spec <-  ugarchspec(variance.model=list(garchOrder=c(1,1)),
                                mean.model = list(armaOrder=arma.order), distribution.model="std")
fit.list$garch11_std <- ugarchfit(spec=garch11_std.spec, data=s)

garch22_std.spec <-  ugarchspec(variance.model=list(garchOrder=c(2,2)),
                                mean.model = list(armaOrder=arma.order), distribution.model="std")

fit.list$garch22_std <- ugarchfit(spec=garch22_std.spec, data=s)

igarch11_std.spec <-  ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                                 mean.model = list(armaOrder=arma.order, include.mean=TRUE),
                                 distribution.model="std")

fit.list$igarch11_std <- ugarchfit(spec=igarch11_std.spec, data=s)

garch_m_11_std.spec  <-  ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                                    mean.model = list(armaOrder=arma.order,
                                                      include.mean=TRUE, archm=TRUE, archpow = 2),
                                    distribution.model="std")
fit.list$garch_m_11_std <- ugarchfit(spec=garch_m_11_std.spec, data=s)

aparch11_std.spec = ugarchspec(variance.model = list(model="apARCH", garchOrder=c(1,1)),
                                mean.model = list(armaOrder=arma.order, include.mean=TRUE),
                                distribution.model="std")

fit.list$aparch11_std <- ugarchfit(spec=aparch11_std.spec, data=s)

info.mat <- sapply(fit.list, infocriteria)
rownames(info.mat) <-  rownames(infocriteria(fit.list[[1]]))

info.mat

for (model in fit.list) {
  variance_model <- model@model$modeldesc$vmodel
  type_dist <- model@model$modeldesc$distribution
  print(paste(variance_model, type_dist))
  print(round(model@fit$matcoef, 4))
}

fit.list$arch1@fit$matcoef
fit.list$garch11@fit$matcoef
round(fit.list$garch11_std@fit$matcoef, 4)

garch11_var <- fit.list$garch11@fit$var   # save the estimated conditional variances
garch11_res2 <- (fit.list$garch11@fit$residuals)^2   # save the estimated squared residuals
plot(garch11_res2, type = "l")
lines(garch11_var, col = "green")

plot(fit.list$garch11_std, which=8) #8:   Empirical Density of Standardized Residuals
plot(fit.list$garch11, which=8) #8:   Empirical Density of Standardized Residuals
plot(fit.list$garch11, which=9) #9:   QQ-Plot of Standardized Residuals
plot(fit.list$garch11_std, which=9) #9:   QQ-Plot of Standardized Residuals

garchforecast1 <- ugarchforecast(fit.list$garch11_std, data= s, n.ahead = 1, )
plot(garchforecast1, which=1) #1:   Time Series Prediction (unconditional)
plot(garchforecast1, which=3) #3:   Sigma Prediction (unconditional)

# ---- GARCH ----
modelspec = ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                       mean.model = list(armaOrder=c(0,0)),
                       distribution.model="norm")

fit = ugarchfit(modelspec, s)
show(fit)
summary(fit)

# ar1 is the AR1 coefficient of the mean model (here very small and basically insignificant),
# alpha1 is the coefficient to the squared residuals in the GARCH equation and
# beta1 is the coefficient to the lagged variance.

# ---- iGARCH ----
modelspec <-  ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(2,2)),
                         mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                         distribution.model="norm")

fit <-  ugarchfit(modelspec, s)
show(fit)
summary(fit)

# ---- GARCH-M ----
modelspec  <-  ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                       mean.model = list(armaOrder=c(0,0), include.mean=TRUE, archm=TRUE, archpow = 2),
                       distribution.model="norm")

fit  <-  ugarchfit(modelspec, log_returnEDP)
show(fit)

# ---- APARCH ----

#aparch_model <- garchFit(formula = ~ aparch(1,1), log_returnGALP, include.delta=TRUE, cond.dist = "std")

#The "signbias" Sign Bias Test of Engle and Ng, if is needed APARCH
# Weighted Ljung-Box and ARCH-LM statistics

modelspec = ugarchspec(variance.model = list(model="apARCH", garchOrder=c(1,1)),
                       mean.model = list(armaOrder=c(0,0), include.mean=TRUE, archm=TRUE, archpow = 2),
                       distribution.model="norm")

fit = ugarchfit(modelspec, log_returnEDP)
show(fit)

# Let's plot the "squared residuals" and the "estimated conditional variance":
ug_var <- ugfit@fit$var   # save the estimated conditional variances
ug_res2 <- (ugfit@fit$residuals)^2   # save the estimated squared residuals
plot(ug_res2, type = "l")
lines(ug_var, col = "green")


fwdCast = ugarchforecast(fta,   n.ahead=5,  n.roll=5)
fwdCast
plot(fwdCast, which="all")


# Now let’s use the rugarch’s standard functionality to use this estimated model
# to produce rolling forecasts of $\sigma_t$ and plot them versus $∣log_return∣$
spec           <- getspec(garchfit)
setfixed(spec) <- as.list(coef(garchfit))
garchforecast1 <- ugarchforecast(spec, n.ahead = 1, n.roll = 2499, data = SPYRet_xts["2007-02-01/"], out.sample = 2500)
plot(garchforecast1, which = 4)
# It looks like to me this model does a pretty good job of sensing how long a volatility
# spike will remain elevated, or rather modeling the path of a volatility spike back down
# to long-run mean levels. Since all econometric models use past values to predict current
# values, it cannot foresee the initial spike up in volatility.



