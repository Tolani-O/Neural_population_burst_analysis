functions {
  matrix getPartials(matrix Omega, int reg_index) {
    int ind = rows(Omega);
    matrix[ind, ind] pcor = diag_matrix(rep_vector(1, ind));
    pcor[reg_index,reg_index] = 0;
    for (x in 2:ind) {
      if (x!=reg_index) {
        for (y in 1:(x-1)) {
          if (y!=reg_index) {
            real val = (Omega[x,y]-(Omega[x,reg_index]*Omega[y,reg_index]))/
            (sqrt(1-square(Omega[x,reg_index]))*sqrt(1-square(Omega[y,reg_index])));
            pcor[x,y] = val;
            pcor[y,x] = val;
          }
        }
      }
    }
    return(pcor);
  }
}
data {
  int<lower=1> Nctr;
  int<lower=1> Ncr;
  int<lower=1> c;
  int<lower=1> t;
  int<lower=1> r;
  int<lower=1, upper=r> reg_index;
  vector[Nctr] X_ctr_data;
  vector[Ncr] X_cr_data;
  int<lower=1> X_ctr_info[Nctr,4];
  vector[Nctr] X_ctr_sd;
  vector[Ncr] X_cr_sd;
}
transformed data {
  vector[Nctr] beta_ctr = X_ctr_sd;
  vector[Ncr] sigma_cr = X_cr_sd;
}
parameters {
  vector[c*t*r] Z_ctr_raw;
  vector[Ncr] Z_cr_raw;
  cholesky_factor_corr[r] L_Omega; //Cholesky factorization of Area Correlation
  vector<lower=0>[c*r] gamma_cr_raw; // for P_ct_structured
}
transformed parameters {
  vector[Nctr] P_ctr_unstructured;
  vector[Ncr] P_cr_unstructured;
  matrix[c, r] gamma_c = to_matrix(gamma_cr_raw, c, r);
  {
    matrix[t, r] Z_prime;
    int conf = 0;
    matrix[c*t, r] Z_ctr = to_matrix(Z_ctr_raw, c*t, r);
    P_cr_unstructured = X_cr_data + sigma_cr .* Z_cr_raw;
    for (i in 1:Nctr) {
      int info[4] = X_ctr_info[i];
      int this_c = info[1];
      int this_t = info[2];
      int this_r = info[3];
      if (this_c!=conf) {
        Z_prime = (diag_pre_multiply(gamma_c[this_c], L_Omega)*Z_ctr[(this_c*t-t+1):(this_c*t)]')';
        conf = conf+1;
      }
      P_ctr_unstructured[i] = P_cr_unstructured[info[4]] + Z_prime[this_t, this_r];
    }
  }
}
model {
  Z_ctr_raw ~ std_normal();
  Z_cr_raw ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(1);
  gamma_cr_raw ~ cauchy(0,10);
  X_cr_data ~  normal(P_cr_unstructured, sigma_cr);
  X_ctr_data ~ normal(P_ctr_unstructured, beta_ctr);
}
generated quantities {
 corr_matrix[r] Omega = multiply_lower_tri_self_transpose(L_Omega);
 matrix[r,r] Pcorr= getPartials(Omega, reg_index);
}
