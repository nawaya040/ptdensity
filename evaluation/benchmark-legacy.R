# Running time, peak memory and held-out log density of ptdensity and of the
# 2020 SMCMP software (smc.apt) on the same data and settings.
# One configuration per R process:
#   Rscript evaluation/benchmark-legacy.R <tool> <d> <n> <particles> <max_depth> <seed> <out.csv> [legacy_lib]
# tool is "ptdensity" or "legacy". legacy_lib is a library holding SMCMP built
# from https://github.com/nawaya040/SMCMP (commit 558806e).
args <- commandArgs(TRUE)
tool <- args[1]; d <- as.integer(args[2]); n <- as.integer(args[3])
M <- as.integer(args[4]); depth <- as.integer(args[5]); seed <- as.integer(args[6]); out <- args[7]
.libPaths(c(if (length(args) >= 8) args[8], "C:/Users/naway/project/ptdensity/.local-library", .libPaths()))

# Data on the unit cube: two-component mixture of independent beta variables.
simulate_data <- function(n, d) {
  z <- runif(n) < .5
  x <- matrix(0, n, d)
  for (j in seq_len(d)) x[, j] <- ifelse(z, rbeta(n, 4, 12), rbeta(n, 12, 4))
  x
}
set.seed(seed)
train <- simulate_data(n, d); test <- simulate_data(5000, d)
peak_ws <- function() {
  cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())
  as.numeric(system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE)) / 2^20
}
base_mb <- peak_ws()
set.seed(seed + 1)
t0 <- proc.time()[[3]]
if (tool == "ptdensity") {
  library(ptdensity)
  fit <- density_pt(train, support = cbind(rep(0, d), rep(1, d)),
                    prior = density_prior(max_depth = depth, split_grid = 31, min_node = 5, eta = .1),
                    control = density_control(particles = M, ess_threshold = .1))
  dens <- predict(fit, test)
} else {
  library(SMCMP)
  fit <- smc.apt(train, test, max.resol = depth, n.particle = M, n.grid.L = 31,
                 eta.R = .1, min.obs = 5)
  dens <- fit$pred.density
}
elapsed <- proc.time()[[3]] - t0
row <- data.frame(tool = tool, d = d, n = n, particles = M, max_depth = depth, seed = seed,
                  seconds = elapsed, peak_mb = peak_ws(), start_mb = base_mb,
                  mean_log_density = mean(log(dens)))
write.table(row, out, sep = ",", append = file.exists(out), col.names = !file.exists(out),
            row.names = FALSE)
print(row)
