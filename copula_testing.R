setwd("brms")
devtools::load_all()

bdata_bin = data.frame(n = c(2,1,4,5,1), d = c(5,5,5,5,5),
                       y = rnorm(n = 5,mean = 4,sd=10), c = c(2,1,3,2,4))

mn <- bf(y ~ 1, family = gaussian())
mp <- bf(c ~ 1, family = poisson())
mb <- bf(n | trials(d) ~ 1, family = binomial())

mod <- mn + mp + mb + set_rescor(rescor = TRUE, copula = "gaussian")

t <- brm(mod, data = bdata_bin)
