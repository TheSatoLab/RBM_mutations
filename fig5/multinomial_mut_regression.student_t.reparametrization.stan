

data {
  int K;
  int N;
  int D;
  real S_w;
  matrix[N,2] TS;
  matrix[K,D] COEF;
  int Y[N,K];
  real generation_time;
  int bin_size;
  int Y_sum[N];
}

transformed data {
  real Max_date;
  matrix[N,2] TS_norm;
  Max_date = TS[N,2];
  TS_norm = TS / Max_date;
}

parameters {
  vector[K-1] b0_raw;
  vector[K] b1_n;
  real<lower=0> b1_g;
  // vector[K] b1_raw;
  vector[D] w;
  real<lower=0> S_b;
}

transformed parameters {
  vector[1] Zero;
  vector[K] b0;
  vector[K] b1_raw;
  vector[K] b1;
  vector[K] X;
  matrix[K,2] b;
  matrix[K,N] mu;
  Zero[1] = 0;
  b0 = append_row(Zero,b0_raw);
  X = COEF * w;
  b1_raw = b1_n / sqrt(b1_g);
  b1 = (X + b1_raw ) * S_b;
  b = append_col(b0,b1);
  mu = b * TS_norm';
}


model {
  S_b ~ student_t(4, 0, 1);
  w ~ double_exponential(0,S_w);
  // b1_raw ~ student_t(4, X, 1);
  b1_g ~ gamma(2.0, 2.0);
  b1_n ~ std_normal();

  for (n in 1:N)
    Y[n] ~ multinomial_logit(mu[,n]);

}



generated quantities {
  vector[K] growth_rate;
  vector[D] growth_gain;
  int Y_predict[N,K];
  
  growth_rate = exp(((b1 / Max_date) / bin_size)  * generation_time);
  growth_gain = exp(((w / Max_date) / bin_size)  * generation_time);

  for(n in 1:N){
    Y_predict[n] = multinomial_rng(softmax(mu[,n]),Y_sum[n]);
  }

}






