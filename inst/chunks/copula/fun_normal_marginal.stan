  /**
  * Normal marginal
  *
  * Standardized normal marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * See https://github.com/spinkney/helpful_stan_functions for
  * further documentation
  *
  * @copyright Ethan Alt, Sean Pinkney 2021 \n
  *
  * @param y matrix of normal outcomes
  * @param mu_glm row_vector of regression means
  * @param matrix vector of outcome SD's
  * @return 2D array of matrices containing the random variables
  *         and jacobian adjustments
  */
  matrix[] normal_marginal(matrix y, matrix mu_glm, vector sigma) {
    int N = rows(mu_glm);
    int J = cols(mu_glm);
    array[2] matrix[N, J] rtn;
    // Initialise the jacobian adjustments to zero, as vectorised lpdf will be used
    rtn[2] = rep_matrix(0, N, J);

    for (j in 1 : J) {
      rtn[1][ : , j] = (y[ : , j] - mu_glm[ : , j]) / sigma[j];
      rtn[2][1, j] = normal_lpdf(y[ : , j] | mu_glm[ : , j], sigma[j]);
    }

    return rtn;
  }
