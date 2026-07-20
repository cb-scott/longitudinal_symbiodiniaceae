// Bernoulli (binary) version of coral_dirichlet_partition_site.stan
// Response Y is binary: 0 = no switching, 1 = any switching event
// Likelihood: Bernoulli logit (replaces cumulative ordinal logit)
// Random effects: STUDY/SITE/COLONY_ID nesting (3 levels), parallel to the
// cumulative-outcome site model -- this replaces the old 2-level
// STUDY / STUDY:COLONY_ID structure to match:
//   shifted ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus +
//     (1 | STUDY/SITE/COLONY_ID) +
//     (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec,       cov = A)) +
//     (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec_trait, cov = T_mat)) +
//     (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid)
//
// All other structure unchanged: Dirichlet partition prior on species-level
// variance across phylogenetic, residual, and trait components.
// Key change: sd_1 (phylogenetic/clean.spec), sd_2 (residual/clean.spec_resid),
// and sd_3 (trait/clean.spec_trait) are no longer sampled independently.
// Instead, for each coefficient k (1=Intercept, 2=PC1, 3=PC2, 4=PC3):
//   - total_var[k] ~ prior on total species-level variance
//   - pi_k ~ Dirichlet(1,1,1) flat prior on partition (phylo, resid, trait)
//   - sd_1[k] = sqrt(total_var[k] * pi_k[1])  <- derived, not sampled
//   - sd_2[k] = sqrt(total_var[k] * pi_k[2])  <- derived, not sampled
//   - sd_3[k] = sqrt(total_var[k] * pi_k[3])  <- derived, not sampled
// sd_4 (STUDY), sd_5 (STUDY:SITE), and sd_6 (STUDY:SITE:COLONY_ID) are
// unchanged -- intercept-only, not part of the species-level variance decomposition.

functions {
  /* compute correlated group-level effects with known covariance matrix */
  matrix scale_r_cor_cov(matrix z, vector SD, matrix L, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    matrix[rows(L), cols(L)] LC = diag_pre_multiply(SD, L);
    int rows_z = rows(z);
    int rows_L = rows(L);
    for (icov in 1:rows(Lcov)) {
      for (jcov in 1:icov) {
        if (Lcov[icov, jcov] > 1e-10) {
          for (i in 1:rows_L) {
            for (j in 1:i) {
              int k = (rows_L * (icov - 1)) + i;
              int l = (rows_L * (jcov - 1)) + j;
              r[k] = r[k] + Lcov[icov, jcov] * LC[i, j] * z_flat[l];
            }
          }
        }
      }
    }
    return to_matrix(r, cols(z), rows(z), 0);
  }

  /* compute correlated group-level effects (no covariance matrix) */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> Y;  // binary: 0 = no switch, 1 = any switch
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=1> Kc;
  // group 1: clean.spec (phylogenetic, has covariance matrix A)
  int<lower=1> N_1;
  int<lower=1> M_1;
  array[N] int<lower=1> J_1;
  matrix[N_1, N_1] Lcov_1;
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  vector[N] Z_1_3;
  vector[N] Z_1_4;
  int<lower=1> NC_1;
  // group 2: clean.spec_resid (no covariance matrix)
  int<lower=1> N_2;
  int<lower=1> M_2;
  array[N] int<lower=1> J_2;
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  vector[N] Z_2_3;
  vector[N] Z_2_4;
  int<lower=1> NC_2;
  // group 3: clean.spec_trait (has covariance matrix T_mat)
  int<lower=1> N_3;
  int<lower=1> M_3;
  array[N] int<lower=1> J_3;
  matrix[N_3, N_3] Lcov_3;
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  vector[N] Z_3_3;
  vector[N] Z_3_4;
  int<lower=1> NC_3;
  // group 4: STUDY (intercept only)
  int<lower=1> N_4;
  int<lower=1> M_4;
  array[N] int<lower=1> J_4;
  vector[N] Z_4_1;
  // group 5: STUDY:SITE (intercept only)
  int<lower=1> N_5;
  int<lower=1> M_5;
  array[N] int<lower=1> J_5;
  vector[N] Z_5_1;
  // group 6: STUDY:SITE:COLONY_ID (intercept only)
  int<lower=1> N_6;
  int<lower=1> M_6;
  array[N] int<lower=1> J_6;
  vector[N] Z_6_1;
  int prior_only;
}

transformed data {
  matrix[N, Kc] Xc;
  vector[Kc] means_X;
  for (i in 1:Kc) {
    means_X[i] = mean(X[, i + 1]);  // skip column 1 (intercept handled separately)
    Xc[, i] = X[, i + 1] - means_X[i];
  }
}

parameters {
  vector[Kc] b;
  real Intercept;  // single intercept on log-odds scale (replaces ordered thresholds)

  // ---- NEW: direct Dirichlet partition parameterization ----
  // total_var[k] = total species-level variance for coefficient k
  // (k=1: Intercept, k=2: PC1, k=3: PC2, k=4: PC3)
  vector<lower=0>[M_1] total_var;

  // pi_k = partition of total_var[k] into (phylo, resid, trait)
  // sampled directly as simplexes -- this IS the prior on the partition
  simplex[3] pi_1;  // Intercept partition
  simplex[3] pi_2;  // PC1 partition
  simplex[3] pi_3;  // PC2 partition
  simplex[3] pi_4;  // PC3 partition

  // z and L parameters unchanged -- only the SD scaling changes
  matrix[M_1, N_1] z_1;
  cholesky_factor_corr[M_1] L_1;
  matrix[M_2, N_2] z_2;
  cholesky_factor_corr[M_2] L_2;
  matrix[M_3, N_3] z_3;
  cholesky_factor_corr[M_3] L_3;

  // STUDY, STUDY:SITE, and STUDY:SITE:COLONY_ID SDs unchanged
  vector<lower=0>[M_4] sd_4;
  array[M_4] vector[N_4] z_4;
  vector<lower=0>[M_5] sd_5;
  array[M_5] vector[N_5] z_5;
  vector<lower=0>[M_6] sd_6;
  array[M_6] vector[N_6] z_6;
}

transformed parameters {
  // ---- DERIVED: sd_1, sd_2, sd_3 computed from total_var and partition ----
  // sd_1 = phylogenetic (clean.spec)
  // sd_2 = residual    (clean.spec_resid)
  // sd_3 = trait       (clean.spec_trait)
  vector<lower=0>[M_1] sd_1;
  vector<lower=0>[M_2] sd_2;
  vector<lower=0>[M_3] sd_3;

  {
    // build a matrix of partitions: rows = coefficients, cols = (phylo, resid, trait)
    matrix[M_1, 3] pi_mat;
    pi_mat[1] = pi_1';
    pi_mat[2] = pi_2';
    pi_mat[3] = pi_3';
    pi_mat[4] = pi_4';

    for (k in 1:M_1) {
      sd_1[k] = sqrt(total_var[k] * pi_mat[k, 1]);  // phylogenetic
      sd_2[k] = sqrt(total_var[k] * pi_mat[k, 2]);  // residual
      sd_3[k] = sqrt(total_var[k] * pi_mat[k, 3]);  // trait
    }
  }

  // actual group-level effects -- identical to original brms code
  matrix[N_1, M_1] r_1;
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  vector[N_1] r_1_3;
  vector[N_1] r_1_4;

  matrix[N_2, M_2] r_2;
  vector[N_2] r_2_1;
  vector[N_2] r_2_2;
  vector[N_2] r_2_3;
  vector[N_2] r_2_4;

  matrix[N_3, M_3] r_3;
  vector[N_3] r_3_1;
  vector[N_3] r_3_2;
  vector[N_3] r_3_3;
  vector[N_3] r_3_4;

  vector[N_4] r_4_1;
  vector[N_5] r_5_1;
  vector[N_6] r_6_1;

  // prior contributions
  real lprior = 0;

  // compute actual group-level effects using derived SDs
  r_1 = scale_r_cor_cov(z_1, sd_1, L_1, Lcov_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  r_1_3 = r_1[, 3];
  r_1_4 = r_1[, 4];

  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_1 = r_2[, 1];
  r_2_2 = r_2[, 2];
  r_2_3 = r_2[, 3];
  r_2_4 = r_2[, 4];

  r_3 = scale_r_cor_cov(z_3, sd_3, L_3, Lcov_3);
  r_3_1 = r_3[, 1];
  r_3_2 = r_3[, 2];
  r_3_3 = r_3[, 3];
  r_3_4 = r_3[, 4];

  r_4_1 = (sd_4[1] * (z_4[1]));
  r_5_1 = (sd_5[1] * (z_5[1]));
  r_6_1 = (sd_6[1] * (z_6[1]));

  // ---- PRIORS (lprior) ----
  // Single intercept: normal(-3, 1) reflects low baseline switching probability
  // on log-odds scale (logistic(-3) ≈ 0.05, i.e. ~5% baseline switching rate)
  lprior += normal_lpdf(Intercept | -3, 1);

  // Fixed effects
  lprior += normal_lpdf(b | 0, 1);

  // total_var: prior on the total species-level variance per coefficient
  // half-normal(0,1) on the SD scale = student_t(3,0,1) truncated at 0
  // using lognormal to match your earlier choice and allow heavier tail
  lprior += lognormal_lpdf(total_var | log(1), 1);

  // STUDY, STUDY:SITE, and STUDY:SITE:COLONY_ID SDs -- unchanged from original
  lprior += student_t_lpdf(sd_4 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_5 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_6 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);

  // LKJ priors on correlation matrices -- unchanged
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 1);
  lprior += lkj_corr_cholesky_lpdf(L_3 | 1);
}

model {
  // ---- DIRECT DIRICHLET PRIORS ON PARTITION ----
  // Dirichlet(1,1,1) = flat/uniform over the simplex
  // genuinely uninformative about which source dominates
  // no implicit penalty on asymmetric partitions
  pi_1 ~ dirichlet(rep_vector(1.0, 3));
  pi_2 ~ dirichlet(rep_vector(1.0, 3));
  pi_3 ~ dirichlet(rep_vector(1.0, 3));
  pi_4 ~ dirichlet(rep_vector(1.0, 3));

  // likelihood
  if (!prior_only) {
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;  // single intercept
    mu += Xc * b;
    for (n in 1:N) {
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n]
             + r_1_3[J_1[n]] * Z_1_3[n] + r_1_4[J_1[n]] * Z_1_4[n]
             + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n]
             + r_2_3[J_2[n]] * Z_2_3[n] + r_2_4[J_2[n]] * Z_2_4[n]
             + r_3_1[J_3[n]] * Z_3_1[n] + r_3_2[J_3[n]] * Z_3_2[n]
             + r_3_3[J_3[n]] * Z_3_3[n] + r_3_4[J_3[n]] * Z_3_4[n]
             + r_4_1[J_4[n]] * Z_4_1[n]
             + r_5_1[J_5[n]] * Z_5_1[n]
             + r_6_1[J_6[n]] * Z_6_1[n];
    }
    target += bernoulli_logit_lpmf(Y | mu);  // Bernoulli logit likelihood
  }

  // priors
  target += lprior;

  // z priors unchanged
  target += std_normal_lpdf(to_vector(z_1));
  target += std_normal_lpdf(to_vector(z_2));
  target += std_normal_lpdf(to_vector(z_3));
  target += std_normal_lpdf(z_4[1]);
  target += std_normal_lpdf(z_5[1]);
  target += std_normal_lpdf(z_6[1]);
}

generated quantities {
  // intercept on original (mean-centered predictor) scale
  real b_Intercept = Intercept + dot_product(means_X, b);

  // correlation matrices
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;

  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }

  // ---- VARIANCE PARTITION: directly exposed as generated quantities ----
  // h2[k, 1] = phylogenetic fraction for coefficient k
  // h2[k, 2] = residual fraction for coefficient k
  // h2[k, 3] = trait fraction for coefficient k
  // k: 1=Intercept, 2=PC1, 3=PC2, 4=PC3
  matrix[M_1, 3] h2;
  h2[1, 1] = pi_1[1]; h2[1, 2] = pi_1[2]; h2[1, 3] = pi_1[3];
  h2[2, 1] = pi_2[1]; h2[2, 2] = pi_2[2]; h2[2, 3] = pi_2[3];
  h2[3, 1] = pi_3[1]; h2[3, 2] = pi_3[2]; h2[3, 3] = pi_3[3];
  h2[4, 1] = pi_4[1]; h2[4, 2] = pi_4[2]; h2[4, 3] = pi_4[3];

  // also expose derived SDs for diagnostics/comparison with original model
  vector[M_1] sd_phylo = sd_1;
  vector[M_2] sd_resid = sd_2;
  vector[M_3] sd_trait = sd_3;
}
