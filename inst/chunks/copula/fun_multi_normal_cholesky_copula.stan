  /**
  * Multivariate Normal Copula log density (Cholesky parameterisation)
  *
  * See https://github.com/spinkney/helpful_stan_functions for
  * further documentation
  *
  * @copyright Ethan Alt, Sean Pinkney 2021
  *
  * @param U Matrix of outcomes from marginal calculations
  * @param L Cholesky factor of the correlation/covariance matrix
  * @return log density of distribution
  */
  real multi_normal_cholesky_copula_lpdf(matrix U, matrix L) {
    int N = rows(U);
    int J = cols(U);
    matrix[J, J] Gammainv = chol2inv(L);
    return -N * sum(log(diagonal(L))) - 0.5 * sum(add_diag(Gammainv, -1.0) .* crossprod(U));
  }
