setwd("brms")
devtools::load_all()

n_out = mnormt::rmnorm(n = 500, mean = c(0, 0),
                       varcov = cbind(c(1,-0.6),c(-0.6,1)))
bdata = data.frame(
  y1 = n_out[,1],
  y2 = n_out[,2],
  c = rpois(500, 10)
)

bdata_bin = data.frame(n = c(2,1,4,5,1), d = c(5,5,5,5,5),
                       y1 = n_out[,1],
                       y2 = n_out[,2],
                       c = c(2,1,3,2,4),
                       c2 = c(10,5,1,8,20))


mn <- bf(y1 ~ 1, family = gaussian())
mn2 <- bf(y2 ~ 1, family = gaussian())
mp <- bf(c ~ 1, family = poisson())
mp2 <- bf(c2 ~ 1, family = poisson())
mb <- bf(n | trials(d) ~ 1, family = binomial())

mod_n <- mn + mp + mn2 + set_rescor(rescor = TRUE, copula = "gaussian")

t_n500 <- brm(mod_n, data = bdata, cores = 4, backend = "cmdstanr")
