  /**
  * Mixed Copula Log-Probability Function
  *
  * @copyright Andrew Johnson 2022 \n
  *
  * @param marginals Nested arrays of matrices from marginal calculations
  * @param L Cholesky Factor Correlation
  * @return Real log-probability
  */
  real centered_gaussian_copula_cholesky_(matrix[,] marginals, matrix L, int[] order) {
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
    return multi_normal_cholesky_copula_lpdf(U[ : , order] | L) + adj;
  }
