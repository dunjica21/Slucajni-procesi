
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

bs_call_scalar <- function(S, K, T, r, sigma) {
  if (sigma <= 0) return(max(0, S - K * exp(-r*T)))
  d1 <- (log(S/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S * pnorm(d1) - K * exp(-r*T) * pnorm(d2)
}


mjd_call <- function(S, K, T, r, sigma, lambda, m, delta, N = 50L) {
  k_tilde <- exp(m + 0.5 * delta^2) - 1
  n_vals  <- 0:(N - 1L)
  
 
  pois <- dpois(n_vals, lambda = lambda * T)
  
 
  sum(vapply(n_vals, function(n) {
    sigma_n <- sqrt(sigma^2 + (n / T) * delta^2)
    S_n     <- S * exp(-lambda * k_tilde * T + (m + 0.5 * delta^2) * n)
    bs_call_scalar(S_n, K, T, r, sigma_n)
  }, numeric(1)) * pois)
}


impvol <- function(price, S, K, T, r,
                   lower = 1e-8, upper = 5, tol = 1e-8, maxit = 100) {
  f <- function(sig) bs_call_scalar(S, K, T, r, sig) - price
  
  intrinsic <- max(0, S - K * exp(-r*T))
  upper_bd  <- S
  if (price < intrinsic - 1e-12 || price > upper_bd + 1e-12) return(NA_real_)
  
 
  u <- upper; fu <- f(u); l <- lower; fl <- f(l)
  tries <- 0
  while (sign(fu) == sign(fl) && tries < 15) {
    u <- u * 1.5
    fu <- f(u); tries <- tries + 1
  }
  if (sign(fu) == sign(fl)) return(NA_real_)
  
  uniroot(f, c(l, u), tol = tol, maxiter = maxit)$root
}


smile_df <- function(Ks, S, T, r, sigma, lambda, m, delta, N = 50L) {
  prices <- vapply(Ks, function(K)
    mjd_call(S, K, T, r, sigma, lambda, m, delta, N = N), numeric(1))
  ivs <- mapply(impvol, price = prices, K = Ks,
                MoreArgs = list(S = S, T = T, r = r))
  tibble(K = Ks, iv = ivs)
}


S0     <- 100
Tmat   <- 1.0
r      <- 0.05
sigma0 <- 0.20
lambda0<- 0.10
m0     <- -0.10
delta0 <- 0.30

Ks <- seq(60, 140, by = 2)  


df1 <- smile_df(Ks, S0, Tmat, r, sigma0, lambda0, m0, delta0)
p1 <- ggplot(df1, aes(K, iv)) +
  geom_line(size = 1.05) +
  geom_hline(yintercept = sigma0, linetype = "dashed") +
  labs(title = "Volatilni osmeh u Mertonovom modelu",
       x = " K", y = "Implicirani volatilitet (Black–Scholes)")

print(p1)


sigmas <- c(0.20, 0.22, 0.24)
df2 <- bind_rows(lapply(sigmas, function(sg) {
  smile_df(Ks, S0, Tmat, r, sg, lambda0, m0, delta0) %>%
    mutate(serija = paste0("sigma=", sprintf("%.0f%%", 100*sg)))
}))
p2 <- ggplot(df2, aes(K, iv, color = serija)) +
  geom_line(size = 1.05) +
  labs(title = "Uticaj osnovne volatilnosti na oblik osmeha u Mertonovom modelu",
       
       x = " K", y = "Implicirani volatilitet", color = NULL)

print(p2)


m_vals <- c(-0.1, 0.0, 0.1)
df3 <- bind_rows(lapply(m_vals, function(mv) {
  smile_df(Ks, S0, Tmat, r, sigma0, lambda0, mv, delta0) %>%
    mutate(serija = paste0("m=", sprintf("%.1f", mv)))
}))
p3 <- ggplot(df3, aes(K, iv, color = serija)) +
  geom_line(size = 1.05) +
  labs(title = "Uticaj prosečne veličine skoka (m) na MJD osmeh",
       subtitle = sprintf("S0=%g, T=%.1f, r=%.0f%%, sigma=%.0f%%, lambda=%.0f%%, delta=%.1f",
                          S0, Tmat, 100*r, 100*sigma0, 100*lambda0, delta0),
       x = " K", y = "Implicirani volatilitet", color = NULL)

print(p3)


deltas <- c(0.3, 0.5, 0.7)
df4 <- bind_rows(lapply(deltas, function(dlt) {
  smile_df(Ks, S0, Tmat, r, sigma0, lambda0, m0, dlt) %>%
    mutate(serija = paste0("delta=", sprintf("%.1f", dlt)))
}))
p4 <- ggplot(df4, aes(K, iv, color = serija)) +
  geom_line(size = 1.05) +
  labs(title = "Uticaj volatilnosti skoka (delta) na MJD osmeh",
       subtitle = sprintf("S0=%g, T=%.1f, r=%.0f%%, sigma=%.0f%%, lambda=%.0f%%, m=%.2f",
                          S0, Tmat, 100*r, 100*sigma0, 100*lambda0, m0),
       x = " K", y = "Implicirani volatilitet", color = NULL)

print(p4)

