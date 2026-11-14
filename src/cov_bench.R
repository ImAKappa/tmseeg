library(CovTools)
library(bench)

# 1. Benchmark execution speed and memory allocation
timing <- bench::mark(
  AIRM = CovMean(covs, method = "AIRM"),
  LERM = CovMean(covs, method = "LERM"),
  check = FALSE, # AIRM and LERM produce slightly different matrices
  iterations = 20
)
print(timing[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

# 2. Compare numerical differences between the two mean matrices
mean_airm <- CovMean(covs, method = "AIRM")
mean_lerm <- CovMean(covs, method = "LERM")

# Relative Frobenius norm difference
rel_frob_diff <- norm(mean_airm - mean_lerm, type = "F") / norm(mean_airm, type = "F")

# Geodesic (AIRM) distance between the two estimates
geo_dist <- CovDist(array(c(mean_airm, mean_lerm), dim = c(61, 61, 2)), method = "AIRM")[1, 2]

cat(sprintf("Relative Frobenius difference: %.6f\n", rel_frob_diff))
cat(sprintf("Riemannian (AIRM) distance:    %.6f\n", geo_dist))


min_eigs <- apply(covs, 3, function(S) min(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
summary(min_eigs)
