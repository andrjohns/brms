
  real inv_Phi_log_lambda(real log_p) {
    vector[8] a
        = [3.3871328727963666080e+00, 1.3314166789178437745e+02,
           1.9715909503065514427e+03, 1.3731693765509461125e+04,
           4.5921953931549871457e+04, 6.7265770927008700853e+04,
           3.3430575583588128105e+04, 2.5090809287301226727e+03]';
    vector[7] b
        = [4.2313330701600911252e+01, 6.8718700749205790830e+02,
           5.3941960214247511077e+03, 2.1213794301586595867e+04,
           3.9307895800092710610e+04, 2.8729085735721942674e+04,
           5.2264952788528545610e+03]';
    vector[8] c
        = [1.42343711074968357734e+00, 4.63033784615654529590e+00,
           5.76949722146069140550e+00, 3.64784832476320460504e+00,
           1.27045825245236838258e+00, 2.41780725177450611770e-01,
           2.27238449892691845833e-02, 7.74545014278341407640e-04]';
    vector[7] d
        = [2.05319162663775882187e+00, 1.67638483018380384940e+00,
           6.89767334985100004550e-01, 1.48103976427480074590e-01,
           1.51986665636164571966e-02, 5.47593808499534494600e-04,
           1.05075007164441684324e-09]';
    vector[8] e
        = [6.65790464350110377720e+00, 5.46378491116411436990e+00,
           1.78482653991729133580e+00, 2.96560571828504891230e-01,
           2.65321895265761230930e-02, 1.24266094738807843860e-03,
           2.71155556874348757815e-05, 2.01033439929228813265e-07]';
    vector[7] f
        = [5.99832206555887937690e-01, 1.36929880922735805310e-01,
           1.48753612908506148525e-02, 7.86869131145613259100e-04,
           1.84631831751005468180e-05, 1.42151175831644588870e-07,
           2.04426310338993978564e-15]';
    real log_q = log_p <= log(0.5) ? log_diff_exp(log(1), log_sum_exp(log_p, log(0.5))) : log_diff_exp(log_p, log(0.5));
    int log_q_sign = log_p <= log(0.5) ? -1 : 1;
    real val;

    if (log_p == log(1)) {
      return positive_infinity();
    }

    if (log_q <= log(.425)) {
      real log_r;
      real log_rtn;
      real log_agg_a;
      real log_agg_b;
      vector[8] log_a = log(a);
      vector[7] log_b = log(b);
      log_r = log_diff_exp(log(.180625), 2 * log_q);
      log_agg_a = log_sum_exp(log_a[8] + log_r, log_a[7]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[6]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[5]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[4]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[3]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[2]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[1]);

      log_agg_b = log_sum_exp(log_b[7] + log_r, log_b[6]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[5]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[4]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[3]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[2]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[1]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, 0);

      log_rtn = log_q + log_agg_a - log_agg_b;
      return log_q_sign * exp(log_rtn);
    } else {
      real log_r = log_q_sign == -1 ? log_p : log_diff_exp(log(1), log_p);

      if (is_inf(log_r)) {
        return 0;
      }

      log_r = log(sqrt(-log_r));

      if (log_r <= log(5.0)) {
        vector[8] log_c = log(c);
        vector[7] log_d = log(d);
        real log_agg_c;
        real log_agg_d;
        log_r = log_diff_exp(log_r, log(1.6));

        log_agg_c = log_sum_exp(log_c[8] + log_r, log_c[7]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[6]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[5]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[4]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[3]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[2]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[1]);

        log_agg_d = log_sum_exp(log_d[7] + log_r, log_d[6]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[5]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[4]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[3]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[2]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[1]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, 0);

        val = exp(log_agg_c - log_agg_d);
      } else {
        vector[8] log_e = log(e);
        vector[7] log_f = log(f);
        real log_agg_e;
        real log_agg_f;
        log_r = log_diff_exp(log_r, log(5));

        log_agg_e = log_sum_exp(log_e[8] + log_r, log_e[7]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[6]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[5]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[4]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[3]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[2]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[1]);

        log_agg_f = log_sum_exp(log_f[7] + log_r, log_f[6]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[5]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[4]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[3]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[2]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[1]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, 0);

        val = exp(log_agg_e - log_agg_f);
      }
      if (log_q_sign == -1)
        return -val;
    }
    return val;
  }

  real inv_Phi_log_fun(real log_p) {
    real log_BIGINT = log(2000000000);
    return log_p >= log(0.9999) ? -inv_Phi_log_lambda(
               log_diff_exp(log_BIGINT, log_BIGINT + log_p) - log_BIGINT)
                       : inv_Phi_log_lambda(log_p);
  }
