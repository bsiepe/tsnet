## code to prepare `ts_data` dataset goes here
set.seed(2023)

# Random parameters
mod1 <- graphicalVAR::randomGVARmodel(6, probKappaEdge = .4, probBetaEdge = .3)
mod2 <- graphicalVAR::randomGVARmodel(6, probKappaEdge = .5, probBetaEdge = .2)

# Simulate data
dat1 <- graphicalVAR::graphicalVARsim(250, beta = mod1$beta, kappa = mod1$kappa)
dat2 <- graphicalVAR::graphicalVARsim(250, beta = mod2$beta, kappa = mod2$kappa)

# Scale
dat1 <- as.data.frame(apply(dat1, 2, scale))
dat2 <- as.data.frame(apply(dat2, 2, scale))

# Add indicator
dat1$id <- rep("ID1", 250)
dat2$id <- rep("ID2", 250)

# Combine
ts_data <- as.data.frame(rbind(dat1, dat2))

usethis::use_data(ts_data, overwrite = TRUE)
