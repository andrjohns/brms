setwd("brms")
devtools::load_all()

bdata_bin = data.frame(n = c(2,1,4,5,1), d = c(5,5,5,5,5),
                       y = rep(1, 5), c = rep(1,5))

mn <- bf(y ~ 1, family = gaussian())
mp <- bf(c ~ 1, family = poisson())
mb <- bf(n | trials(d) ~ 1, family = binomial())

mod <- mn + mp + mb + set_rescor(rescor = TRUE, copula = "gaussian")

make_stancode(mod, data = bdata_bin)
