  /**
  * Bernoulli marginal
  *
  * Bernoulli marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * See https://github.com/spinkney/helpful_stan_functions for
  * further documentation
  *
  * @copyright Ethan Alt, Sean Pinkney 2021 \n
  *
  * @param y int[,] 2d array of binary outcomes
  * @param mu_glm matrix of regression means
  * @param u_raw matrix of nuisance latent variables
  * @return 2D array of matrices containing the random variables
  *         and jacobian adjustments
  */
  matrix[] bernoulli_marginal(int[,] y, matrix mu_glm, matrix u_raw) {
    int N = rows(mu_glm);
    int J = cols(mu_glm);
    matrix[N, J] matrix_y = to_matrix(y);
    matrix[N, J] mu_glm_logit = 1 - inv_logit(mu_glm);

    matrix[N, J] Lbound = matrix_y .* mu_glm_logit;
    matrix[N, J] UmL = fabs(matrix_y - mu_glm_logit);

    return {inv_Phi(Lbound + UmL .* u_raw), log(UmL)};
  }
