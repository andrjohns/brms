  real std_normal_log_qf2_scal(real log_p) {
    real log_a[8]
        = {1.2199838032983212, 4.8914137334471356, 7.5865960847956080,
            9.5274618535358388, 10.734698580862359, 11.116406781896242,
            10.417226196842595, 7.8276718012189362};
    real log_b[8] = {0.,
                                    3.7451021830139207,
                                    6.5326064640478618,
                                    8.5930788436817044,
                                    9.9624069236663077,
                                    10.579180688621286,
                                    10.265665328832871,
                                    8.5614962136628454};
    real log_c[8]
        = {0.3530744474482423, 1.5326298343683388, 1.7525849400614634,
            1.2941374937060454, 0.2393776640901312, -1.419724057885092,
            -3.784340465764968, -7.163234779359426};
    real log_d[8] = {0.,
                                    0.7193954734947205,
                                    0.5166395879845317,
                                    -0.371400933927844,
                                    -1.909840708457214,
                                    -4.186547581055928,
                                    -7.509976771225415,
                                    -20.67376157385924};
    real log_e[8]
        = {1.8958048169567149, 1.6981417567726154, 0.5793212339927351,
            -1.215503791936417, -3.629396584023968, -6.690500273261249,
            -10.51540298415323, -15.41979457491781};
    real log_f[8] = {0.,
                                    -0.511105318617135,
                                    -1.988286302259815,
                                    -4.208049039384857,
                                    -7.147448611626374,
                                    -10.89973190740069,
                                    -15.76637472711685,
                                    -33.82373901099482};

    real val;
    real LOG_HALF = -log2();
    real log_q = log_p <= LOG_HALF ? log_diff_exp(LOG_HALF, log_p)
                                      : log_diff_exp(log_p, LOG_HALF);
    int log_q_sign = log_p <= LOG_HALF ? -1 : 1;

    if (log_p == negative_infinity()) {
      return negative_infinity();
    }
    if (log_p == 0) {
      return positive_infinity();
    }

    if (log_q <= -0.85566611005772) {
      real log_r = log_diff_exp(-1.71133222011544, 2 * log_q);
      real log_agg_a = log_sum_exp(log_a[8] + log_r, log_a[7]);
      real log_agg_b = log_sum_exp(log_b[8] + log_r, log_b[7]);

      for (i in 0:5) {
        log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[6 - i]);
        log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[6 - i]);
      }

      return log_q_sign * exp(log_q + log_agg_a - log_agg_b);
    } else {
      real log_r = log_q_sign == -1 ? log_p : log1m_exp(log_p);

      if (is_inf(log_r)) {
        return 0;
      }

      log_r = log(sqrt(-log_r));

      if (log_r <= 1.60943791243410) {
        real log_agg_c;
        real log_agg_d;
        log_r = log_diff_exp(log_r, 0.47000362924573);
        log_agg_c = log_sum_exp(log_c[8] + log_r, log_c[7]);
        log_agg_d = log_sum_exp(log_d[8] + log_r, log_d[7]);

        for (i in 0:5) {
          log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[6 - i]);
          log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[6 - i]);
        }
        val = exp(log_agg_c - log_agg_d);
      } else {
        real log_agg_e;
        real log_agg_f;
        log_r = log_diff_exp(log_r, 1.60943791243410);
        log_agg_e = log_sum_exp(log_e[8] + log_r, log_e[7]);
        log_agg_f = log_sum_exp(log_f[8] + log_r, log_f[7]);

        for (i in 0:5) {
          log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[6 - i]);
          log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[6 - i]);
        }
        val = exp(log_agg_e - log_agg_f);
      }
      if (log_q_sign == -1)
        return -val;
    }
    return val;
  }

  matrix std_normal_log_qf2(matrix log_p_mat) {
    int R = rows(log_p_mat);
    int C = cols(log_p_mat);
    matrix[R, C] rtn;

    for (c in 1:C) {
      for (r in 1:R) {
        rtn[r, c] = std_normal_log_qf2_scal(log_p_mat[r, c]);
      }
    }

    return rtn;
  }
