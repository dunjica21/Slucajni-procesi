a     <- 0.8
b     <- 0.04
sigma <- 0.12
t0    <- 0

rho <- sqrt(a^2 + 2 * sigma^2)

cir_phi <- function(theta, a, sigma) {
  rho <- sqrt(a^2 + 2 * sigma^2)
  num <- 2 * rho * exp((a + rho) * theta / 2)
  den <- (a + rho) * (exp(rho * theta) - 1) + 2 * rho
  num / den
}

cir_g <- function(theta, a, sigma) {
  rho <- sqrt(a^2 + 2 * sigma^2)
  2 * (exp(rho * theta) - 1) / ((a + rho) * (exp(rho * theta) - 1) + 2 * rho)
}

cir_yield_theorem <- function(t, T, rt, a, b, sigma) {
  theta <- pmax(T - t, 0)
  y <- numeric(length(theta))
  nz <- theta > 0
  if (any(nz)) {
    ph <- cir_phi(theta[nz], a, sigma)
    gg <- cir_g(theta[nz], a, sigma)
    y[nz] <- (gg / theta[nz]) * rt - (2 * a * b / sigma^2) * (log(ph) / theta[nz])
  }
  y[!nz] <- rt
  y
}

T <- seq(1e-4, 100, by = 0.1)
r0_vals <- c(0.015, 0.035, 0.065)

yields <- lapply(r0_vals, function(r0) cir_yield_theorem(t0, T, r0, a, b, sigma))

R_inf <- 2 * a * b / (a + rho)

boje <- c("red", "blue", "darkgreen")

ylim <- range(unlist(yields), R_inf)
plot(T, yields[[3]], type = "l", lwd = 2, ylim = ylim, col = boje[3],
     xlab = "T", ylab = "R(0,T)", main = "CIR krive prinosa")
lines(T, yields[[2]], lwd = 2, col = boje[2])
lines(T, yields[[1]], lwd = 2, col = boje[1])
abline(h = R_inf, lty = 3)

legend("topright",
       legend = c(paste0("r0 = ", r0_vals[3]),
                  paste0("r0 = ", r0_vals[2]),
                  paste0("r0 = ", r0_vals[1]),
                  paste0("R\u221E ≈ ", sprintf("%.3f", R_inf))),
       lwd = c(2,2,2,1), lty = c(1,1,1,3), col = c(boje[3], boje[2], boje[1], "black"),
       bty = "n")
