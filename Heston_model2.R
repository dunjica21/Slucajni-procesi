heston_cf <- function(s0, v0, vbar, kappa, sigma_v, r, rho, t, w) {
  i  <- 1i
  alpha <- -0.5 * w * w - 0.5 * i * w
  beta  <- kappa - rho * sigma_v * i * w
  gamma <- 0.5 * sigma_v^2
  h     <- sqrt(beta * beta - 4 * alpha * gamma)
  rplus  <- (beta + h) / sigma_v^2
  rminus <- (beta - h) / sigma_v^2
  g      <- rminus / rplus
  C <- kappa * (rminus * t - (2 / sigma_v^2) * log((1 - g * exp(-h * t)) / (1 - g)))
  D <- rminus * (1 - exp(-h * t)) / (1 - g * exp(-h * t))
  exp(C * vbar + D * v0 + i * w * log(s0 * exp(r * t)))
}

heston_call_cf <- function(s0, v0, vbar, kappa, sigma_v, r, rho, t, K,
                           upper = 100, rel.tol = 1e-8) {
  i <- 1i
  # π1 integral
  integrand1 <- function(w) {
    w <- as.numeric(w)
    Re( exp(-i * w * log(K)) *
          heston_cf(s0, v0, vbar, kappa, sigma_v, r, rho, t, w - i) /
          (i * w * heston_cf(s0, v0, vbar, kappa, sigma_v, r, rho, t, -i)) )
  }
  pi1 <- 0.5 + (1/pi) * integrate(integrand1, lower = 0, upper = upper,
                                  rel.tol = rel.tol)$value
  
  integrand2 <- function(w) {
    w <- as.numeric(w)
    Re( exp(-i * w * log(K)) *
          heston_cf(s0, v0, vbar, kappa, sigma_v, r, rho, t, w) / (i * w) )
  }
  pi2 <- 0.5 + (1/pi) * integrate(integrand2, lower = 0, upper = upper,
                                  rel.tol = rel.tol)$value
  
  s0 * pi1 - exp(-r * t) * K * pi2
}

heston_call_cf_vec <- function(s0, v0, vbar, kappa, sigma_v, r, rho, t, K,
                               upper = 100, rel.tol = 1e-8) {
  sapply(K, function(k) heston_call_cf(s0, v0, vbar, kappa, sigma_v, r, rho, t, k,
                                       upper = upper, rel.tol = rel.tol))
}

ex_val <- heston_call_cf(s0 = 1, v0 = 0.16, vbar = 0.16,
                         kappa = 1, sigma_v = 2, r = 0, rho = -0.8,
                         t = 10, K = 2)
ex_val

