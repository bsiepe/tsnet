## Code to create example fit object
data(ts_data)
ts_data1 <- ts_data[1:100,1:3]
ts_data2 <- ts_data[250:349,1:3]

fit1 <- tsnet::stan_gvar(data = ts_data1,
                        n_chains = 2,
                        n_cores = 1)
fit2 <- tsnet::stan_gvar(data = ts_data2,
                        n_chains = 2,
                        n_cores = 1)

# Obtain samples
samples1 <- tsnet::stan_fit_convert(fit1,
                                    return_params = c("beta", "pcor"))
samples2 <- tsnet::stan_fit_convert(fit2,
                                    return_params = c("beta", "pcor"))

fit_data <- list(samples1, samples2)
usethis::use_data(fit_data, overwrite = TRUE)
