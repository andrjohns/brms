  /**
  * Binomial marginal
  *
  * Binomial marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * See https://github.com/spinkney/helpful_stan_functions for
  * further documentation
  *
  * @copyright Ethan Alt, Sean Pinkney 2021 \n
  *
  * @param num int[,] 2D array of numerator integers
  * @param den int[,] 2D array of denominator integers
  * @param mu_glm matrix of regression means
  * @param u_raw matrix of nuisance latent variables
  * @return 2D array of matrices containing the random variables
  *         and jacobian adjustments
  */
  matrix[] binomial_marginal(int[,] num, int[,] den, matrix mu_glm,
                                  matrix u_raw) {
    int N = rows(mu_glm);
    int J = cols(mu_glm);
    matrix[N, J] mu_glm_logit = inv_logit(mu_glm);
    matrix[N, J] rtn[2];

    for (j in 1 : J) {
      for (n in 1 : N) {
        real Ubound = binomial_cdf(num[n, j], den[n, j], mu_glm_logit[n, j]);
        real Lbound = 0;
        real UmL;
        if (num[n, j] > 0) {
          Lbound = binomial_cdf(num[n, j] - 1, den[n, j], mu_glm_logit[n, j]);
        }
        UmL = Ubound - Lbound;
        rtn[1][n, j] = inv_Phi(Lbound + UmL * u_raw[n, j]);
        rtn[2][n, j] = log(UmL);
      }
    }

    return rtn;
  }
