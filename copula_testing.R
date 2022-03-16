setwd("~/brms")
devtools::load_all()

set.seed(2022)

cor <- round(randcorr::randcorr(4), 1) # Rounding for clarity

# Generate correlated standard-normal variates
n_out <- mnormt::rmnorm(n = 500, mean = c(0, 0, 0, 0),
                        varcov = cor)

# Put on uniform [0, 1] scale using standard normal CDF (Phi function)
p_out <- pnorm(n_out)

# Generate denominators for binomial
d <- sample(10:100, 500, replace = T)

# Transform uniform variates to desired distributions using
# quantile function with desired parameters
bdata <- data.frame(
  y1 = qnorm(p_out[,1], 10, 4),
  y2 = qnorm(p_out[,2], 3, 6),
  c = qpois(p_out[,3], 15),
  n = qbinom(p_out[,4], prob = 0.3, size = d),
  d = d
)

colnames(cor) <- c("y1","y2","c","n")
rownames(cor) <- c("y1","y2","c","n")


mn <- bf(y1 ~ 1, family = gaussian())
mp <- bf(c ~ 1, family = poisson())
mn2 <- bf(y2 ~ 1, family = gaussian())
mb <- bf(n | trials(d) ~ 1, family = binomial())
mod_n <- mn + mp + mn2 + mb + set_rescor(rescor = TRUE, copula = "gaussian")

t_n <- brm(mod_n, data = bdata, cores = 4)
