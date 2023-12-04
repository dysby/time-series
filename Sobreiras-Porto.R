# Imports ----
library(tidyverse, quietly = TRUE)
library(fpp3, quietly = TRUE)
library(imputeTS)
library(latex2exp)
library(patchwork)
library(ggpubr)

# Load ----
qualityAR03 <- readxl::read_excel("./data/QualidadeARO3.xlsx") %>%
  ts(.,start=c(2020,1,1),frequency=24*366) %>%
  as_tsibble() %>%
  rename(Location = key, Ozone = value) %>%
  mutate(index = as_datetime(round_date(index, unit = "hour")))

y <- qualityAR03 %>% filter(Location=="Sobreiras-Porto")

# GGplot ----
y %>% gg_tsdisplay(Ozone, plot_type = "partial", lag_max = 200)

# The Sobreira series show variance not constant and some breaking behavior in June, therefore a Box-Cox transformation will be applied to stabilize the series.

# Transform ----
lambda <- y %>%
  features(Ozone, features = guerrero) %>%
  pull(lambda_guerrero)

lambda

#lambda = 0 # log transform

y %>%
  gg_tsdisplay(box_cox(Ozone, lambda), plot_type = "partial", lag_max = 200) +
  labs(y = "",
       title = TeX(paste0(
         "Difference h=24 Box-Cox Transformation $\\lambda$ = ",
         round(lambda,2))))

# The Sobreira series present an ACF slowly decaying because of the presence of a seasonal component at lags 24 (daily seasonal), which motivates to apply a seassonal diferencing.

y %>%
  gg_tsdisplay(difference(box_cox(Ozone, lambda), lag = 24), plot_type = "partial", lag_max = 200) +
  #gg_tsdisplay(box_cox(Ozone + 1, lambda=lambda), plot_type = "partial", lag_max = 50) +
  labs(y = "",
       title = TeX(paste0(
         "Difference h=24 Box-Cox Transformation $\\lambda$ = ",
         #round(lambda,2), ", $\\log(Ozone + 1)$")))
         round(lambda,2))))

# The seasonal difference series present a pattern resembling random noise with
# constant variance, although there are some points with large values around
# september. ACF presents a pattern of exponential decay along
# the lags, but still have a peak at lag 24 meaning that the seasonal component
# still have a role. Since the PACf shows a exponential decaying pattern in
# seasonal lags, the seasonal component present have the characteristics of a
# SMA(1) (or even SMA(2) if we consider the peak at the lag 48 in the ACF being relevant).

# Stationary tests ----

adf <- function(x, ...) {
  out <- tseries::adf.test(x, k = 1)
  c(adf_stat = unname(out$statistic), adf_pvalue = out$p.value)
}

tseries::adf.test(box_cox(y$Ozone, lambda), k=1)

y %>%
  mutate(Ozone = box_cox(Ozone, lambda)) %>%
  features(Ozone, list(unitroot_kpss, adf, unitroot_ndiffs, ~ unitroot_nsdiffs(., .period = 24)))

# difference and remove 1st NA value
tseries::adf.test(difference(box_cox(y$Ozone, lambda), lag = 24)[-c(1:24)], k=1)
tseries::adf.test(difference(box_cox(y$Ozone, lambda))[-1], k=1)

# Diff ACF plots ----
# The acf plot of (d=1, D=1) have only one peak at lag 24, the PACF have decreasing behaviour at the seasonal lags but perhaps should have a faster  be faster decreasing pattern. The seasonal component of the time series is (0, 1, 1)[24].

p7 <- y %>% ACF(Ozone, lag_max = 200) %>%
  autoplot() + labs(title = "Original")
p8 <- y %>% PACF(Ozone, lag_max = 200) %>% autoplot()

p5 <- y %>% mutate(diff = difference(Ozone)) %>%
  ACF(diff, lag_max = 200) %>% autoplot() + labs(title = "Difference series")
p6 <- y %>% mutate(diff = difference(Ozone)) %>%
  PACF(diff, lag_max = 200) %>% autoplot()

p1 <- y %>% mutate(diff_S = difference(Ozone, lag=24)) %>%
  ACF(diff_S, lag_max = 200) %>% autoplot() +
  labs(title = "Seasonal (24) difference series")
p2 <- y %>% mutate(diff_S = difference(Ozone, lag=24)) %>%
  PACF(diff_S, lag_max = 200) %>% autoplot()

p3 <- y %>% mutate(diff_SN = difference(difference(Ozone, lag=24))) %>%
  ACF(diff_SN, lag_max = 200) %>% autoplot() + labs(title = "Seasonal (24) and\n non-seasonal differences series")
p4 <- y %>% mutate(diff_SN = difference(difference(Ozone, lag=24))) %>%
  PACF(diff_SN, lag_max = 200) %>% autoplot()

# (p7 + p8) / (p5 + p6) / (p1 + p2) / (p3 + p4) #/
(p7 + p8) / (p1 + p2)

y %>% mutate(diff_SN = difference(box_cox(Ozone, lambda), lag=24)) %>% gg_tsdisplay(diff_SN, plot_type = "partial", lag_max = 200)
# plot TS  diff 24
y %>% mutate(diff_SN = difference(difference(box_cox(Ozone, lambda), lag=24))) %>% gg_tsdisplay(diff_SN, plot_type = "partial", lag_max = 200)

y %>% mutate(diff_SN = difference(difference(box_cox(Ozone, lambda), lag=24), lag=24)) %>% gg_tsdisplay(diff_SN, plot_type = "partial", lag_max = 200)

# sugere SARIMA (2,0,?)(0, 1, 1)
fit <- y %>%
  model(
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() +
  ggtitle(paste(fit$Location, format(fit %>% pull(2)))) + ylim(-0.1,0.5)
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot() + ylim(-0.1,0.5)
p1 + p2

# the SARIMA (2,0,0)(0, 1, 1) model show mostly random white noise pattern in ACF and PACF of the residuals.
# The residuals do not have seasonal component dependencies. AS seen by the partial and autocovariance at lags 24, 48, 72, etc being all zero.
# We don't need to test an alternative having one more SMA term to account for the lag 48 peak in the ACF of seasonal difference series.

fit <- y %>%
  model(
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 1, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() +
  ggtitle(paste(fit$Location, format(fit %>% pull(2)))) + ylim(-0.1,0.5)
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot() + ylim(-0.1,0.5)
p1 + p2

# The ACF and PACF of the residuals of the SARIMA (2,0,0)(0, 1, 2) model, are not
# similar to white noise that the previous model candidate.
# The fit using SARIMA (2,0,0)(0, 1, 1) and some variations, to evaluate measures AIC, AICc and BIC.

# compare s101 with s012

fit <- y %>%
  model(
    sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima100011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima201011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima300011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima200012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 2, period=24)),
  )

m <- fit %>% select("Location", all_of("sarima200011"))
m <- fit %>% select("Location", all_of("sarima201011"))
#m <- fit %>% select("Location", all_of("sarima300011"))
p1 <- m %>% augment() %>% ACF(.innov, lag_max = 100) %>% autoplot() +
  ggtitle(paste(m$Location, format(m %>% pull(2)))) +
  ylim(-0.1,0.5)
p2 <- m %>% augment() %>% PACF(.innov, lag_max = 100) %>% autoplot() +
  ylim(-0.1,0.5)
p1 + p2

# Model statistics ----

fit %>% pivot_longer(!Location, names_to = "Model name",
                     values_to = "Orders") %>%
  mutate("Model" = format(Orders)) %>%
  glance() %>%
  select(Location, Model, AIC, AICc, BIC)

# The model with the best overall metrics is the one predicted by examining the ACF and PACF (<ARIMA(2,0,0)(0,1,1)[24]>).


# bind_rows( fit %>% glance(),
           # fit2 %>% glance())

# Coefficients ----

fit %>% select(sarima200011) %>% coef()

# The coefficients of the model are all statisticaly relevants (not zero).


# residuals ----
fit %>% gg_tsresiduals(lag_max = 200)

# bind_rows( fit %>% glance(),
#            fit2 %>% glance())

# The Ljung Box test, along the lags until 10, results in p-values all above the decision threshold 0.05 (except the ones at lag 2 and 6, but close to the threshold) indicates that we can consider the residuals uncorrelated.

sapply(c(1:10), function(x) {m %>% features(.innov, ljung_box, lag=x) %>% pull(lb_pvalue)}) %>%
  as_tibble() %>%
ggplot(aes(x=1:10, y=value)) +
  geom_bar(stat='identity', width =0.05) +
  geom_point() + ylim(0,1) + ylab("p-value") +
  geom_abline(intercept = .05, slope = 0, lty = 3) +
  scale_x_discrete(limits=1:10) +
  xlab("h") +
  ggtitle("p-value of the Ljung-Box test")

# Create QQplot and Fitter vs Residuals with ggplot2 package

p1 <- fit %>% select(sarima200011) %>% augment() %>%
  mutate(.innov = scale(.innov)) %>%
  ggqqplot(".innov",
         facet.by = ".model",
         title = "Standardized Residuals QQ Plot")

p2 <- fit %>% select(sarima200011) %>% augment() %>%
  ggplot(aes(x = .fitted, y = .resid)) +
  geom_point() +
  labs(title = "Fitted vs Residuals")
p1 + p2

# The normal QQ plot fo the standardized residuals show deviations from normality
# for extreme values. The residuals do not appear to have a Normal distribution,
# and to use the model for forcasting we must be were of this limitation and use
# bootstrapping when forecasting to get random values of innovations with the same distribution as the residuals.
# The Residuals - Fitted plot shows a scatted pattern without clear indication
# of correlations or increase of variance with higher values.

# Model Selection ----

# Save the city object
saveRDS(fit, file = "models/sobreiras-porto_fit_s011.rds")
# Load the city object as city
fit <- readRDS("models/sobreiras-porto_fit_s011.rds")

best_model_name = "sarima200011"
#best_model_name = "sarima400012"

# Forecast ----

knitr::kable(fit %>% select("Location", all_of(best_model_name)) %>%
               forecast(h=5) %>%
               hilo(level = c(95)),
             caption = "Forecast 5 time periods ahead (95\\% CI)")

fit %>% select("Location", all_of(best_model_name)) %>%
  forecast(h=5, bootstrap=TRUE) %>%
  autoplot(y %>%
             filter(between(index,as_datetime("2020-12-15"), as_datetime("2021-01-01"))),
           level = 95,
           alpha = 0.5)

## Model Performance ----

fit %>% select("Location", all_of(best_model_name)) %>% accuracy()

### Cross Validation ----

# https://otexts.com/fpp3/training-test.html
training <- y %>% filter(month(index) <= 11)
test <- y %>% filter(month(index) > 11)

train_fit <- training %>%
  model(
    sarimafinal = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24))
  )
#train_fit <- fit[best_model_name][[1]][[1]]$data %>% model(sarima = fit[best_model_name][[1]][[1]]$model$formula)

# DONOT RUN - very slow
# fits5 <- fitted(train_fit, h = 5)
# training %>%
#   autoplot(Ozone) +
#   autolayer(fits5, .fitted, col = "#D55E00") +
#   labs(title = "Cross Validation Forecast",
#        y = "Ozone concentrations")

# One-step forecasts on test data
# It is common practice to fit a model using training data,
# and then to evaluate its performance on a test data set
# we now apply the model to the test data
fits5 <- train_fit %>%
  refit(test)

test %>%
  autoplot(Ozone) +
  autolayer(fitted(fits5, h = 1), .fitted, col = "#D55E00") +
  labs(title = "Cross Validation on test data",
       y = "Ozone concentrations",
       x = "Date")

# The Cross Validation plot shows that the model is predicting accurately on test data, having only small errors.

fits5 %>% accuracy()

# The model performance metrics cross validated are equivalent to the measures obtained if we fit with the full series.
