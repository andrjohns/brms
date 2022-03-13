  /* Returns the inverse of the matrix whose Cholesky factor is L
  *
  * Reimplemented for compatibility with Stan < 2.26
  *
  * Args:
  *   L: Matrix that is a Cholesky factor
  * Returns:
  *   The matrix inverse of L * L'
  */
  matrix chol2inv2(matrix L) {
    int K = rows(L);
    matrix[K, K] K_identity = diag_matrix(rep_vector(1.0, K));
    matrix[K, K] L_inv = mdivide_left_tri_low(L, K_identity);
    matrix[K, K] rtn;

    if (K == 0) {
      return L;
    }

    if (K == 1) {
      rtn[1, 1] = inv_square(L[1, 1]);
    } else {
      for (k in 1:K) {
        rtn[k, k] = dot_self(tail(L_inv[ : , k], K + 1 - k));
        for (j in (k + 1):K) {
          int Kmj = K - j;
          rtn[k, j] = dot_product(tail(L_inv[ : , k], Kmj + 1),
                                  tail(L_inv[ : , j], Kmj + 1));
          rtn[j, k] = rtn[k, j];
        }
      }
    }

    return rtn;
  }
