  /**
    * Multivariate Normal Copula log density (Cholesky parameterisation)
    *
    * \f{aligned}{
    * c(\mathbf{u}) &= \frac{\partial^d C}{\partial \mathbf{\Phi_1}\cdots \partial \mathbf{\Phi_d}} \\
    *               &= \prod_{i=1}^{n} \lvert \Sigma \rvert^{-1/2} \exp \left\{-\frac{1}{2}
    *                              \begin{pmatrix} \mathbf{\Phi}^{-1}(u_1) \\ \vdots  \\ \mathbf{\Phi}^{-1}(u_d)  \end{pmatrix}^T
    *                              (\Sigma^{-1} - I)
    *                              \begin{pmatrix} \mathbf{\Phi}^{-1}(u_1) \\ \vdots  \\ \mathbf{\Phi}^{-1}(u_d)  \end{pmatrix}
    *                              \right\} \\
    *                 &= \lvert \Sigma \rvert^{-n/2} \exp\left\{ -\frac{1}{2} \text{tr}\left( \left[ \Sigma^{-1} - I \right] \sum_{i=1}^n \Phi^{-1}(u_i) \Phi^{-1}(u_i)' \right) \right\} \hspace{0.5cm} \\
    *                 &= \lvert L \rvert^{-n} \exp\left\{  -\frac{1}{2} \text{tr}\left( \left[ L^{-T}L^{-1}- I \right] Q'Q  \right) \right\} \\
    *                 &= \lvert L \rvert^{-n} \exp\left\{  -\frac{1}{2} \sum_{j=1}^{d}\sum_{i=1}^{n} \left( \left[ L^{-T}L^{-1}- I \right] \odot Q'Q  \right) \right\} \\
    *  \f}
    * where
    * \f[
    *    Q_{n \times d} = \begin{pmatrix} \mathbf{\Phi}^{-1}(u_1) \\ \vdots  \\ \mathbf{\Phi}^{-1}(u_d)  \end{pmatrix}
    *  \f] 
    * and \f$ \odot \f$ is the elementwise multiplication operator. \n
    * Note that \f$\det \Sigma  = |LL^T| = |L |^2 = \big(\prod_{i=1}^d L_{i,i}\big)^2\f$ so \f$\sqrt{\det \Sigma} = \prod_{i=1}^d L_{i,i}\f$.
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

  /**
    * Normal marginal
    *
    * Standardized normal marginal for mixed continuous-discrete Gaussian
    * copula.
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

  /**
  * Bernoulli marginal
  *
  * Bernoulli marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * The lb and ub:
  *    if y == 0, upper bound at inv_Phi( y, mu )
  *    if y == 1, lower bound at inv_Phi( y - 1, mu )
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

  /**
  * Binomial marginal
  *
  * Binomial marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * The lb and ub:
  *      Always upper bound at inv_Phi( y, mu )
  *      If n != 0, lower bound at inv_Phi(n-1, mu)
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
    array[2] matrix[N, J] rtn;
    
    for (j in 1 : J) {
      for (n in 1 : N) {
        real Ubound = binomial_cdf(num[n, j], den[n, j], mu_glm_logit[n, j]);
        real Lbound = 0;
        if (num[n, j] > 0) {
          Lbound = binomial_cdf(num[n, j] - 1, den[n, j], mu_glm_logit[n, j]);
        }
        real UmL = Ubound - Lbound;
        rtn[1][n, j] = inv_Phi(Lbound + UmL * u_raw[n, j]);
        rtn[2][n, j] = log(UmL);
      }
    }
    
    return rtn;
  }

  /**
  * Poisson marginal
  *
  * Poisson marginal for mixed continuous-discrete Gaussian
  * copula.
  *
  * The lower-bound and upper-bound: \n
  *      The upper-bound is always at \code{.cpp} inv_Phi( y, mu ) \endcode \n
  *      If \f$y \ne 0\f$, lower-bound at \code{.cpp} inv_Phi(y-1, mu) \endcode. \n
  *      At \f$y = 0\f$ the lower-bound is \f$ 0 \f$.
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

  /**
  * Mixed Copula Log-Probability Function
  *
  * @copyright Andrew Johnson 2022 \n
  *
  * @param marginals Nested arrays of matrices from marginal calculations
  * @param L Cholesky Factor Correlation
  * @return Real log-probability
  */
  real centered_gaussian_copula_cholesky_(matrix[,] marginals, matrix L) {
    // Extract dimensions of final outcome matrix
    int N = rows(marginals[1][1]);
    int J = rows(L);
    matrix[N, J] U;
    
    // Iterate through marginal arrays, concatenating the outcome matrices by column
    // and aggregating the log-likelihoods (from continuous marginals) and jacobian
    // adjustments (from discrete marginals)
    real adj = 0;
    int pos = 1;
    for (m in 1 : size(marginals)) {
      int curr_cols = cols(marginals[m][1]);
      
      U[ : , pos : (pos + curr_cols - 1)] = marginals[m][1];
      
      adj += sum(marginals[m][2]);
      pos += curr_cols;
    }
    
    // Return the sum of the log-probability for copula outcomes and jacobian adjustments
    return multi_normal_cholesky_copula_lpdf(U | L) + adj;
  }
