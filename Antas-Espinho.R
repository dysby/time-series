# Imports ----
library(tidyverse, quietly = TRUE)
library(fpp3, quietly = TRUE)
library(imputeTS)
library(latex2exp)
library(patchwork)

# Load ----
qualityAR03 <- readxl::read_excel("./data/QualidadeARO3.xlsx") %>%
  ts(.,start=c(2020,1,1),frequency=24*366) %>%
  as_tsibble() %>%
  rename(Location = key, Ozone = value) %>%
  mutate(index = as_datetime(round_date(index, unit = "hour")))

y <- qualityAR03 %>% filter(Location=="Antas-Espinho")

# GGplot ----
y %>% gg_tsdisplay(Ozone, plot_type = "partial", lag_max = 200)

# Transform ----
lambda <- y %>%
  features(Ozone, features = guerrero) %>%
  pull(lambda_guerrero)

lambda

#lambda = 0 # log transform
y %>%
  gg_tsdisplay(difference(box_cox(Ozone, lambda), lag = 24), plot_type = "partial", lag_max = 200) +
  #gg_tsdisplay(box_cox(Ozone + 1, lambda=lambda), plot_type = "partial", lag_max = 50) +
  labs(y = "",
       title = TeX(paste0(
         "Difference h=24 Box-Cox Transformation $\\lambda$ = ",
         #round(lambda,2), ", $\\log(Ozone + 1)$")))
         round(lambda,2))))

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

# sugere SARIMA (?,0,?)(0, 1, 2)
fit <- y %>%
  model(
    sarima000011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 0) + PDQ(0, 1, 1, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
p1 + p2

fit <- y %>%
  model(
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
p1 + p2

fit <- y %>%
  model(
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 2, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
p1 + p2

fit <- y %>%
  model(
    sarima000011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 1, period=24)),
  )
p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
p1 + p2

# sarima200011 residual ACF and PACF have a small peak in lag 2
# will test for adding more 1,2 AR or 1,2 MA terms

fit <- y %>%
  model(
    sarima000011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima000011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima200011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima201011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima202011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 2) + PDQ(0, 1, 1, period=24)),
    sarima300011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 0) + PDQ(0, 1, 1, period=24)),
    sarima301011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 1) + PDQ(0, 1, 1, period=24)),
    sarima400011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(4, 0, 0) + PDQ(0, 1, 1, period=24)),
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

# fit <- y %>%
#   model(
#     sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 2, period=24)),
#   )
# p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot() + ggtitle(format(fit %>% pull(2)))
# p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
# p1 + p2

# fitAuto ----

fit <- y %>%
  model(
    sarima000011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 0) + PDQ(0, 1, 1, period=24)),
    #sarima = ARIMA(box_cox(Ozone, lambda) ~ pdq(0:2,1,0:2) + PDQ(0, 1, 1:2, period=24),
    #               trace = TRUE, approximation = FALSE)
    # sarima = ARIMA(box_cox(Ozone, lambda) ~ pdq(0,1,1) + PDQ(0, 1, 2, period=24),
    #                trace = TRUE, approximation = FALSE)
    # sarima = ARIMA(box_cox(Ozone, lambda) ~ pdq(1,1,1) + PDQ(0, 1, 1, period=24),
    #                trace = TRUE, approximation = FALSE)
    # sarima200210 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0:3,0,0:3) + PDQ(0:2, 1, 0:1, period=24),
    #                      trace = TRUE, approximation = FALSE, selection_metric = )
    # sarima002012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 0) + PDQ(0, 1, 0, period=24)),
    # sarima200210 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2,0,0) + PDQ(2, 1, 0, period=24)),
    # sarima013010 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0,1,3) + PDQ(0, 1, 0, period=24)),
    #sarima_1__1_ = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0:3,1,0:3) + PDQ(0, 1, 1, period=24),
    #                      trace = TRUE, approximation = FALSE),
    # sarima_1__2_ = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0:2,0,0:2) + PDQ(0:1, 2, 0:1, period=24),
    #                      trace = TRUE, approximation = FALSE),
    # sarima101011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1,0,1) + PDQ(0, 1, 1, period=24),
    #                trace = TRUE, approximation = FALSE)
    # sarima100011 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1,0,0) + PDQ(0, 1, 1, period=24),
    #                      trace = TRUE, approximation = FALSE)
  )

# Best is pdq(2, 1, 1) + PDQ(1, 1, 1, period=24))

# fit <- y %>%
#   model(
#     #arima010011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(0, 1, 0) + PDQ(0, 1, 1, period=24)),
#     #sarima = ARIMA(box_cox(Ozone, lambda),
#     sarima = ARIMA(box_cox(Ozone, lambda) ~ pdq(0:3, 0:1, 0:3) + PDQ(0:1, 1, 0:1, period=24),
#                    trace = TRUE,
#                    approximation = FALSE,
#                    greedy=FALSE,
#                    stepwise=FALSE)
#     #sarima = ARIMA(box_cox(Ozone, lambda) ~ pdq(0:1, 1, 0:1) + PDQ(0:1, 1, 0:1, period=24), trace = TRUE, approximation = FALSE, greedy=FALSE, stepwise=FALSE)
#   )

# fitManual ----

# fit <- y %>%
#   model(
#         #sarima000012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 0) + PDQ(0, 1, 2, period=24)),
#         #sarima100012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 0) + PDQ(0, 1, 2, period=24)),
#         sarima200012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 0) + PDQ(0, 1, 2, period=24)),
#         sarima300012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 0) + PDQ(0, 1, 2, period=24)),
#         #sarima001012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 1) + PDQ(0, 1, 2, period=24)),
#         sarima002012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 2) + PDQ(0, 1, 2, period=24)),
#         sarima003012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(0, 0, 3) + PDQ(0, 1, 2, period=24)),
#         sarima101012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 1) + PDQ(0, 1, 2, period=24)),
#         sarima201012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 1) + PDQ(0, 1, 2, period=24)),
#         sarima301012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 1) + PDQ(0, 1, 2, period=24)),
#         sarima102012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 2) + PDQ(0, 1, 2, period=24)),
#         sarima202012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 2) + PDQ(0, 1, 2, period=24)),
#         sarima302012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 2) + PDQ(0, 1, 2, period=24)),
#         sarima103012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(1, 0, 3) + PDQ(0, 1, 2, period=24)),
#         sarima203012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(2, 0, 3) + PDQ(0, 1, 2, period=24)),
#         sarima303012 = ARIMA(box_cox(Ozone, lambda) ~ 0 + pdq(3, 0, 3) + PDQ(0, 1, 2, period=24)),
#   )
#
# fit2 <- y %>%
#   model(
#     #sarima100011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 0, 0) + PDQ(0, 1, 1, period=24)),
#     #sarima001011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(0, 0, 1) + PDQ(0, 1, 1, period=24)),
#     sarima010011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(0, 1, 0) + PDQ(0, 1, 1, period=24)),
#     sarima110011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 1, 0) + PDQ(0, 1, 1, period=24)),
#     sarima210011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(2, 1, 0) + PDQ(0, 1, 1, period=24)),
#     sarima011011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(0, 1, 1) + PDQ(0, 1, 1, period=24)),
#     sarima111011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 1, 1) + PDQ(0, 1, 1, period=24)),
#     sarima211011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(2, 1, 1) + PDQ(0, 1, 1, period=24)),
#     sarima012011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(0, 1, 2) + PDQ(0, 1, 1, period=24)),
#     sarima112011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 1, 2) + PDQ(0, 1, 1, period=24)),
#     #sarima212011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(2, 1, 2) + PDQ(0, 1, 1, period=24)),
#     #sarima113011 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 1, 3) + PDQ(0, 1, 1, period=24)),
#     #sarima113012 = ARIMA(box_cox(Ozone, lambda) ~ pdq(1, 1, 3) + PDQ(0, 1, 2, period=24)),
#   )

# model(arima000011 = ARIMA(Ozone ~ pdq(0, 1, 0) + PDQ(0, 1, 1, period=24)),
#       arima010011 = ARIMA(Ozone ~ pdq(1, 1, 0) + PDQ(0, 1, 1, period=24)),
#       arima010011 = ARIMA(Ozone ~ pdq(0, 1, 1) + PDQ(0, 1, 1, period=24)),
#       arima101110 = ARIMA(Ozone ~ pdq(1, 0, 1) + PDQ(1, 1, 0, period=24)),
#       arima101011 = ARIMA(Ozone ~ pdq(1, 0, 1) + PDQ(0, 1, 1, period=24)),
#       arima202210 = ARIMA(Ozone ~ pdq(2, 0, 2) + PDQ(2, 1, 0, period=24)))
#gg_tsdisplay(plot_type = "partial", lag_max = 100)

# p1 <- fit %>% augment() %>% ACF(.innov, lag_max = 150) %>% autoplot()
# p2 <- fit %>% augment() %>% PACF(.innov, lag_max = 150) %>% autoplot()
# p1 + p2

# [1] "<ARIMA(2,0,0)(0,1,2)[24]>"

# Model statistics ----

fit %>% pivot_longer(!Location, names_to = "Model name",
                     values_to = "Orders") %>%
  mutate("Model" = format(Orders)) %>%
  glance() %>%
  select(Location, Model, AIC, AICc, BIC)

# bind_rows( fit %>% glance(),
           # fit2 %>% glance())

# Coefficients ----

fit %>% coef()

# residuals ----
fit %>% gg_tsresiduals(lag_max = 200)

bind_rows( fit %>% glance(),
           fit2 %>% glance())

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

library(ggpubr)
ggqqplot(fit %>% augment() %>% mutate(.innov = scale(.innov)),
         ".innov",
         facet.by = ".model",
         title = "Standardized Residuals QQ Plot")

augment(fit) %>%
  ggplot(aes(x = .fitted, y = .resid)) +
  geom_point() +
  labs(x = "Fitted", y = "Residuals")

# Model Selection ----

# Save the city object
saveRDS(fit, file = "models/antas-espinho_fit_s011.rds")
# Load the city object as city
fit <- readRDS("models/antas-espinho_fit_s011.rds")

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

