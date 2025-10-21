a     <- 0.30  
b     <- 0.035  
sigma <- 0.02   

R_inf <- b - sigma^2 / (2 * a^2)

vasicek_yield <- function(T, r0, a, b, sigma) {
  B <- (1 - exp(-a * T)) / a
  A <- exp( (b - sigma^2/(2 * a^2)) * (B - T) - (sigma^2 / (4 * a)) * B^2 )
  y <- - (log(A) - B * r0) / T
  y[T == 0] <- r0
  y
}

T <- seq(0.0001, 100, by = 0.1)
r0_vals <- c(0.015, 0.03, 0.06)


ys <- sapply(r0_vals, function(r0) vasicek_yield(T, r0, a, b, sigma))


boje <- c("red", "blue", "darkgreen")


ylim <- range(c(ys, R_inf))
plot(T, ys[,3], type = "l", lwd = 2, ylim = ylim, col = boje[3],
     xlab = "Vreme T", ylab = "Prinos R(0,T)",
     main = "Vasicek krive prinosa")
lines(T, ys[,2], lwd = 2, col = boje[2])
lines(T, ys[,1], lwd = 2, col = boje[1])
abline(h = R_inf, lty = 3)

legend("topright",
       legend = c(paste0("r0 = ", r0_vals[3]),
                  paste0("r0 = ", r0_vals[2]),
                  paste0("r0 = ", r0_vals[1]),
                  paste0("R∞ = ", sprintf("%.3f", R_inf))),
       lwd = c(2,2,2,1), lty = c(1,1,1,3), col = c(boje[3], boje[2], boje[1], "black"), bty = "n")
