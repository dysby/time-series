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

y <- qualityAR03 %>% filter(Location=="VNTelha-Maia")

# Original Display ----
y %>% gg_tsdisplay(Ozone, plot_type = "partial", lag_max = 200)

# The series present a path with a discontinuity in June and some values
# greater than usual in September.
# It appear to have some non stationary variance, some peaks and breaks, therefore
# a transformation will be applied.
# The ACF shows indication of seasonality at lags 24, 48, 74, etc, motivating a
# differentiation on the seasonal component.
# The series does not show any trend, it appears that there is not need to difference
# in the regular component.

# check components Plot
y %>% filter(index >= "2020-01-01" & index <="2020-01-07" ) %>% model(stl = STL(Ozone ~ season(period = 24))) %>% components() %>% select(season_24) %>% autoplot()

# Using STL decomposition by LOESS, to determine the trend, seasonal and remainder components.
# The trend component does not show an evident increase or decrease pattern on the long term, so the series is stationary in the mean.
# Because the series is very long, the plots do not show clear pattern of seasonality.
# But looking only at a smaller set of points of the series, the STL decomposition plot
# shows a clear pattern corresponding to daily (24 lags) seasonality.
# The series have another seasonal component corresponding to week seasonality.
y %>% filter(index >= "2020-01-01" & index <="2020-01-31" ) %>% model(stl = STL(Ozone)) %>% components() %>% autoplot() + ggtitle("January")
# The SARIMA model can only deal with one seasonal period, the most dominant seasonal
# component is one corresponding to a daily seasonality therefore the SARIMA model
# period will be of 24 lags.


# Transform ----
lambda <- y %>%
  features(Ozone, features = guerrero) %>%
  pull(lambda_guerrero)

lambda

#lambda = 0 # log transform
y %>%
  gg_tsdisplay(difference(box_cox(Ozone, lambda), lag = 24),
               plot_type = "partial", lag_max = 200) +
  labs(y = "",
       title = TeX(paste0(
         "Difference h=24 Box-Cox Transformation $\\lambda$ = ",
         #round(lambda,2), ", $\\log(Ozone + 1)$")))
         round(lambda,2))))

# After seasonal differentiating in the transformed series, it is clear that it
# is stationary and the variance is almost constant, except in the regions the
# original series was inconsistent (September).
# The ACF and PACF show a pattern similar to other Location series, a peak at
# lag 24 (ACF), one small peak at lag 48 (ACF), and decaying peaks ate seasonal lags 24, 48, 72, etc in the PACF.
# Theses patterns indicate a SARIMA model with a SMA(1) order (maybe a SMA(2)
# because the peak at 48 lag in the ACF), and for the regular
# component the first lags are rapidly decaying in the ACF and there are only
# two partial-autocorrelations at lags 1 and 2 with significance.
# It is expected that a SARIMA (2,0,0)(0,1,1)[24] will model well the time
# dependency structure of the series.

# Stationary tests ----

adf <- function(x, ...) {
  out <- tseries::adf.test(x, k = 1)
  c(adf_stat = unname(out$statistic), adf_pvalue = out$p.value)
}

tseries::adf.test(box_cox(y$Ozone, lambda), k=1)

y %>%
  mutate(Ozone = box_cox(Ozone, lambda)) %>%
  features(Ozone, list(unitroot_kpss, adf, unitroot_ndiffs, ~ unitroot_nsdiffs(., .period = 24)))

y %>%
  mutate(Ozone = difference(box_cox(Ozone, lambda))) %>%
  drop_na() %>%
  features(Ozone, list(unitroot_kpss, adf, unitroot_ndiffs, ~ unitroot_nsdiffs(., .period = 24)))

# kpss test (p-value = 0.1) for the seasonal difference transformed series have no evidence to
# reject that it is stationary.

# difference and remove 1st NA value
tseries::adf.test(difference(box_cox(y$Ozone, lambda), lag = 24)[-c(1:24)], k=1)
tseries::adf.test(difference(box_cox(y$Ozone, lambda))[-1], k=1)

# Diff ACF plots ----
# The acf plot of (d=1, D=1) have only one peak at lag 24, the PACF have decreasing behaviour at the seasonal lags but perhaps should have a faster be faster decreasing pattern. The seasonal component of the time series is (0, 1, 1)[24].

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

# The ACF and PACF plots of the residuals of the SARIMA (2,0,0)(0, 1, 1) model are mostly as white noise.
# Having only small correlations around lag 24, corresponding to interactions with
# the regular component. These interactions as considered non significant because
# the values are near the threshold and lags 1, 2, 3, etc are all non significant.
# Adding more AR or AM degrees to the model will not be meaningful or better whiteness similarity of the ACF and PACF.
# The seasonal order selected is correct, the correlations at lags 24, 48, etc (are non significant).

fit <- y %>%
  model(
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima100011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima201011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima300011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 0) + PDQ(0, 1, 1, period=24)),
  )


m <- fit %>% select("Location", all_of("sarima200011"))
m <- fit %>% select("Location", all_of("sarima100011"))
#m <- fit %>% select("Location", all_of("sarima300011"))
p1 <- m %>% augment() %>% ACF(.innov, lag_max = 100) %>% autoplot() +
  ggtitle(paste(m$Location, format(m %>% pull(2)))) +
  ylim(-0.1,0.5)
p2 <- m %>% augment() %>% PACF(.innov, lag_max = 100) %>% autoplot() +
  ylim(-0.1,0.5)
p1 + p2

# fit <- y %>%
#   model(
#     sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 2, period=24)),
#   )
# p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
# p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
# p1 + p2

# Model statistics ----

fit %>% pivot_longer(!Location, names_to = "Model name",
                     values_to = "Orders") %>%
  mutate("Model" = format(Orders)) %>%
  glance() %>%
  select(Location, Model, AIC, AICc, BIC)

# bind_rows( fit %>% glance(),
           # fit2 %>% glance())

# Coefficients ----
fit %>% select(sarima200011) %>% coef()

# The test statistics for the model coefficients have p-values very small/zero
# so they are all significant.

# residuals ----
fit %>% gg_tsresiduals(lag_max = 200)

# bind_rows( fit %>% glance(),
#            fit2 %>% glance())

plot.box.ljung <- function (
    z, k = 10, main = "p-value of the Ljung-Box test", ylab = "p-value"
) {
  p <- rep(NA, k)

  for (i in 1:k) {
    p[i] <- Box.test(z, i, type = "Ljung-Box")$p.value
  }

  ggplot(as_tibble(p), aes(x=1:k, y=p)) +
    geom_bar(stat='identity', width =0.05) +
    geom_point() + ylim(0,1) + ylab(ylab) +
    geom_abline(intercept = .05, slope = 0, lty = 3) +
    ggtitle(main)
  #plot(p, type = 'h', ylim = c(0,1), lwd = 3, main = main, ylab = ylab)
  #abline(h = c(0,.05), lty = 3)
}

for (i in 1:k) {
  p[i] <- Box.test(z, i, type = "Ljung-Box")$p.value
}

plot.box.ljung(fit %>% augment() %>%  pull(.innov), k=10)

# fit %>% augment() %>%  features(.innov, ljung_box, lag=1)
# fit %>% augment() %>%  features(.innov, ljung_box, lag=2)
# fit %>% augment() %>%  features(.innov, ljung_box, lag=3, dof=4)
# fit %>% augment() %>%  features(.innov, ljung_box, lag=4, dof=4)
# fit %>% augment() %>%  features(.innov, ljung_box, lag=5, dof=4)
# fit %>% augment() %>%  features(.innov, ljung_box, lag=6, dof=4)

# Box.test(fit %>% augment() %>% pull(.resid), type = "Ljung-Box", lag = 5)

# qqnorm(fit %>% augment() %>% filter(.model == "sarima102012") %>% pull(.innov) %>% scale())
# qqnorm(fit %>% augment() %>% pull(.innov) %>% scale())
#
# qqnorm(fit2 %>% augment() %>% filter(.model == "sarima011011") %>% pull(.innov))
# qqnorm(fit2 %>% augment() %>% filter(.model == "sarima011011") %>% pull(.innov))

# Create QQplot with ggplot2 package
ggplot(fit %>% augment()) +
  aes(sample = scale(.innov), group=.model, color=.model) +
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle("Normal Q-Q Plot (innov)") +
  xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles")

ggqqplot(fit %>% augment() %>% mutate(.innov = scale(.innov)),
         ".innov",
         facet.by = ".model",
         title = "Standardized Residuals QQ Plot")

augment(fit %>% select(sarima200011)) %>%
  ggplot(aes(x = .fitted, y = .resid)) +
  geom_point() +
  labs(x = "Fitted", y = "Residuals")

# Model Selection ----

# Save the city object
saveRDS(fit, file = "models/vntelha-maia_fit_s012.rds")
# Load the city object as city
fit <- readRDS("models/vntelha-maia_fit_s012.rds")

best_model_name = "sarima201011"
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
    sarimafinal = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 1, period=24))
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
# It is common practice to fit a model using training data, and then to evaluate its performance on a test data set
# we now apply the model to the test data
train_fit %>%
  refit(test) %>%
  accuracy()

# gtsummary::as_kable(gtsummary::tbl_summary(t), format = "latex")


