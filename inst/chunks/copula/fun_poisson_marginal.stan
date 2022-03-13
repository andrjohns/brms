  /**
  * Poisson marginal
  *
  * Poisson marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * See https://github.com/spinkney/helpful_stan_functions for
  * further documentation
  *
  * @copyright Ethan Alt, Sean Pinkney 2021 \n
  *
  * @param y int[,] 2D array of integer counts
  * @param mu_glm matrix of regression means
  * @param u_raw matrix of nuisance latent variables
  * @return 2D array of matrices containing the random variables
  *         and jacobian adjustments
  */
  matrix[] poisson_marginal(int[,] y, matrix mu_glm, matrix u_raw) {
    int N = rows(mu_glm);
    int J = cols(mu_glm);
    matrix[N, J] mu_glm_exp = exp(mu_glm);
    array[2] matrix[N, J] rtn;

    for (j in 1 : J) {
      for (n in 1 : N) {
        real Ubound = poisson_cdf(y[n, j], mu_glm_exp[n, j]);
        real Lbound = 0;
        if (y[n, j] > 0) {
          Lbound = poisson_cdf(y[n, j] - 1, mu_glm_exp[n, j]);
        }
        real UmL = Ubound - Lbound;
        rtn[1][n, j] = inv_Phi(Lbound + UmL * u_raw[n, j]);
        rtn[2][n, j] = log(UmL);
      }
    }

    return rtn;
  }
