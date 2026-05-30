# benchmarks/profile_run.R — profiling workload for SIMD hot-spot identification (G4a)
library(leafblower)
set.seed(42)
n <- 10000L
df <- data.frame(
  age   = sample(c("18-34","35-54","55+"), n, TRUE, c(0.3, 0.4, 0.3)),
  sex   = sample(c("M","F"), n, TRUE),
  edu   = sample(c("HS","Some","BA+"), n, TRUE, c(0.3, 0.3, 0.4)),
  reg   = sample(c("N","S","E","W"), n, TRUE),
  vote  = sample(c("D","R","I"), n, TRUE, c(0.4, 0.35, 0.25))
)
targets <- list(
  age  = c("18-34"=0.30, "35-54"=0.40, "55+"=0.30),
  sex  = c("M"=0.49, "F"=0.51),
  edu  = c("HS"=0.30, "Some"=0.30, "BA+"=0.40),
  reg  = c("N"=0.25, "S"=0.25, "E"=0.25, "W"=0.25),
  vote = c("D"=0.40, "R"=0.35, "I"=0.25)
)
for (i in 1:10) harvest(df, targets, method="oris", max_iterations=500L)
for (i in 1:10) harvest(df, targets, method="raking", max_iterations=500L)
