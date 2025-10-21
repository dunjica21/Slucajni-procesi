d1 <- function(S, K, r, sigma, T) {
  (log(S / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
}

d2 <- function(S, K, r, sigma, T) {
  d1(S, K, r, sigma, T) - sigma * sqrt(T)
}

bs_call <- function(S, K, r, sigma, T) {
  D1 <- d1(S, K, r, sigma, T)
  D2 <- D1 - sigma * sqrt(T)
  S * pnorm(D1) - K * exp(-r * T) * pnorm(D2)
}

bs_put <- function(S, K, r, sigma, T) {
  D1 <- d1(S, K, r, sigma, T)
  D2 <- D1 - sigma * sqrt(T)
  K * exp(-r * T) * pnorm(-D2) - S * pnorm(-D1)
}

delta_call <- function(S, K, r, sigma, T) pnorm(d1(S, K, r, sigma, T))
delta_put  <- function(S, K, r, sigma, T) delta_call(S, K, r, sigma, T) - 1
gamma_bs   <- function(S, K, r, sigma, T) dnorm(d1(S, K, r, sigma, T)) / (S * sigma * sqrt(T))
vega_bs    <- function(S, K, r, sigma, T) S * dnorm(d1(S, K, r, sigma, T)) * sqrt(T)  # per 1.0 sigma
theta_call <- function(S, K, r, sigma, T) {
  D1 <- d1(S, K, r, sigma, T); D2 <- d2(S, K, r, sigma, T)
  -(S * dnorm(D1) * sigma) / (2 * sqrt(T)) - r * K * exp(-r * T) * pnorm(D2)
}
rho_call <- function(S, K, r, sigma, T) K * T * exp(-r * T) * pnorm(d2(S, K, r, sigma, T))

impvol_call <- function(price, S, K, r, T, lower=1e-6, upper=5) {
  f <- function(sg) bs_call(S, K, r, sg, T) - price
  uniroot(f, c(lower, upper))$root
}

K     <- 100
r     <- 0.03
sigma <- 0.20
T     <- 1.0

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

S_grid <- seq(40, 160, length.out = 300)
vols   <- c(0.10, 0.20, 0.40)
df1 <- do.call(rbind, lapply(vols, function(v) {
  data.frame(S=S_grid, C=bs_call(S_grid, K, r, v, T), sigma=factor(v))
}))


p1 <- ggplot(df1, aes(x = S, y = C, color = sigma)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = K, linetype = "dashed", color = "gray50") +
  labs(
    title = "Zavisnost cene kol opcije od cene papira u osnovi (Blek–Šols model)",
    x = "Cena papira u osnovi S",
    y = "Cena kol opcije C(S,t)",
    color = "Volatilnost σ"
  ) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)

