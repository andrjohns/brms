  /* compute monotonic effects
   * Args:
   *   scale: a simplex parameter
   *   i: index to sum over the simplex
   * Returns:
   *   a scalar between 0 and 1
   */
  matrix log_sum_exp_mat(matrix mat1, matrix mat2) {
    int R = rows(mat1);
    int C = cols(mat1);
    matrix[R, C] rtn;

    for (c in 1:C) {
      for (r in 1:R) {
        rtn[r, c] = log_sum_exp(mat1[r, c], mat2[r, c]);
      }
    }

    return rtn;
  }
