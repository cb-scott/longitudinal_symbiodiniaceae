// Ordinal genus-turnover model with explicit host TraitPC x environmental-PC
// interactions. Species-level random intercepts and environmental slopes are
// partitioned between phylogenetic and idiosyncratic components. Trait
// covariance is intentionally omitted because TraitPC1-TraitPC3 enter the
// fixed-effect design matrix directly.

functions {
  matrix scale_r_cor_cov(matrix z, vector SD, matrix L, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    matrix[rows(L), cols(L)] LC = diag_pre_multiply(SD, L);
    int rows_L = rows(L);
    for (icov in 1:rows(Lcov)) {
      for (jcov in 1:icov) {
        if (Lcov[icov, jcov] > 1e-10) {
          for (i in 1:rows_L) {
            for (j in 1:i) {
              int k = (rows_L * (icov - 1)) + i;
              int l = (rows_L * (jcov - 1)) + j;
              r[k] += Lcov[icov, jcov] * LC[i, j] * z_flat[l];
            }
          }
        }
      }
    }
    return to_matrix(r, cols(z), rows(z), 0);
  }

  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}

data {
  int<lower=1> N;
  array[N] int Y;
  int<lower=2> nthres;
  int<lower=1> Kc;
  matrix[N, Kc] X;

  // group 1: host species with phylogenetic covariance
  int<lower=1> N_1;
  int<lower=1> M_1;
  array[N] int<lower=1> J_1;
  matrix[N_1, N_1] Lcov_1;
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  vector[N] Z_1_3;
  vector[N] Z_1_4;
  int<lower=1> NC_1;

  // group 2: residual host-species effects
  int<lower=1> N_2;
  int<lower=1> M_2;
  array[N] int<lower=1> J_2;
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  vector[N] Z_2_3;
  vector[N] Z_2_4;
  int<lower=1> NC_2;

  // group 3: study intercept
  int<lower=1> N_3;
  int<lower=1> M_3;
  array[N] int<lower=1> J_3;
  vector[N] Z_3_1;

  // group 4: site-within-study intercept
  int<lower=1> N_4;
  int<lower=1> M_4;
  array[N] int<lower=1> J_4;
  vector[N] Z_4_1;

  // group 5: colony-within-site intercept
  int<lower=1> N_5;
  int<lower=1> M_5;
  array[N] int<lower=1> J_5;
  vector[N] Z_5_1;

  int prior_only;
}

transformed data {
  matrix[N, Kc] Xc;
  vector[Kc] means_X;
  if (M_1 != 4 || M_2 != 4) {
    reject("Species groups must contain intercept and three environmental slopes.");
  }
  for (i in 1:Kc) {
    means_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - means_X[i];
  }
}

parameters {
  vector[Kc] b;
  ordered[nthres] Intercept;

  // Total unexplained species variance and its phylogenetic/residual split.
  vector<lower=0>[M_1] total_var;
  simplex[2] pi_1;
  simplex[2] pi_2;
  simplex[2] pi_3;
  simplex[2] pi_4;

  matrix[M_1, N_1] z_1;
  cholesky_factor_corr[M_1] L_1;
  matrix[M_2, N_2] z_2;
  cholesky_factor_corr[M_2] L_2;

  vector<lower=0>[M_3] sd_3;
  array[M_3] vector[N_3] z_3;
  vector<lower=0>[M_4] sd_4;
  array[M_4] vector[N_4] z_4;
  vector<lower=0>[M_5] sd_5;
  array[M_5] vector[N_5] z_5;
}

transformed parameters {
  vector<lower=0>[M_1] sd_1;
  vector<lower=0>[M_2] sd_2;
  matrix[N_1, M_1] r_1;
  matrix[N_2, M_2] r_2;
  vector[N_3] r_3_1;
  vector[N_4] r_4_1;
  vector[N_5] r_5_1;
  real lprior = 0;

  {
    matrix[M_1, 2] pi_mat;
    pi_mat[1] = pi_1';
    pi_mat[2] = pi_2';
    pi_mat[3] = pi_3';
    pi_mat[4] = pi_4';
    for (k in 1:M_1) {
      sd_1[k] = sqrt(total_var[k] * pi_mat[k, 1]);
      sd_2[k] = sqrt(total_var[k] * pi_mat[k, 2]);
    }
  }

  r_1 = scale_r_cor_cov(z_1, sd_1, L_1, Lcov_1);
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_3_1 = sd_3[1] * z_3[1];
  r_4_1 = sd_4[1] * z_4[1];
  r_5_1 = sd_5[1] * z_5[1];

  lprior += normal_lpdf(Intercept[1] | 1.81, 1.0);
  lprior += normal_lpdf(Intercept[2] | 3.13, 1.0);
  lprior += normal_lpdf(b | 0, 1);
  lprior += lognormal_lpdf(total_var | log(1), 1);
  lprior += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - M_3 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_4 | 3, 0, 2.5)
    - M_4 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_5 | 3, 0, 2.5)
    - M_5 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 1);
}

model {
  pi_1 ~ dirichlet(rep_vector(1.0, 2));
  pi_2 ~ dirichlet(rep_vector(1.0, 2));
  pi_3 ~ dirichlet(rep_vector(1.0, 2));
  pi_4 ~ dirichlet(rep_vector(1.0, 2));

  if (!prior_only) {
    vector[N] mu = Xc * b;
    for (n in 1:N) {
      mu[n] += r_1[J_1[n], 1] * Z_1_1[n]
             + r_1[J_1[n], 2] * Z_1_2[n]
             + r_1[J_1[n], 3] * Z_1_3[n]
             + r_1[J_1[n], 4] * Z_1_4[n]
             + r_2[J_2[n], 1] * Z_2_1[n]
             + r_2[J_2[n], 2] * Z_2_2[n]
             + r_2[J_2[n], 3] * Z_2_3[n]
             + r_2[J_2[n], 4] * Z_2_4[n]
             + r_3_1[J_3[n]] * Z_3_1[n]
             + r_4_1[J_4[n]] * Z_4_1[n]
             + r_5_1[J_5[n]] * Z_5_1[n];
    }
    for (n in 1:N) {
      target += ordered_logistic_lpmf(Y[n] | mu[n], Intercept);
    }
  }

  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += std_normal_lpdf(to_vector(z_2));
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_4[1]);
  target += std_normal_lpdf(z_5[1]);
}

generated quantities {
  vector[nthres] b_Intercept = Intercept + dot_product(means_X, b);
  matrix[M_1, 2] h2;
  h2[1, 1] = pi_1[1]; h2[1, 2] = pi_1[2];
  h2[2, 1] = pi_2[1]; h2[2, 2] = pi_2[2];
  h2[3, 1] = pi_3[1]; h2[3, 2] = pi_3[2];
  h2[4, 1] = pi_4[1]; h2[4, 2] = pi_4[2];
  vector[M_1] sd_phylo = sd_1;
  vector[M_2] sd_resid = sd_2;
}
