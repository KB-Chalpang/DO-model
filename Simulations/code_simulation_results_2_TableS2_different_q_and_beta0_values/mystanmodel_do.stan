
data {
        int<lower=0> N;                          // number of sites
        array[N] int<lower=0, upper=1> y;        // observed detections
        vector[N] x1;
        vector[N] x2;
        vector[N] x3;
      }
      
      parameters {
        // intercept and ceofficients
        real beta0;
        real beta1;
        real beta2;
        real beta3;
      
        // Visit probability / the composite quantity in DO-M
        real<lower=0.001, upper=1> q;
      
        // Marginal detection probabilities
        real<lower=0, upper=q> theta11;          
        real<lower=0, upper=q> theta10;       
      }
      
      transformed parameters {
        //  detection probabilities
        real p11 = theta11 / q;
        real p10 = theta10 / q;
      
        // occupancy probability
        vector[N] psi = inv_logit(
          beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3
        );
      
        // marginal observation probability
        vector[N] p_obs = theta11 .* psi + theta10 .* (1 - psi);
      }
      
      model {
        // priors
        beta0 ~ normal(0, 2);
        beta1 ~ normal(0, 1);
        beta2 ~ normal(0, 1);
        beta3 ~ normal(0, 1);
        
        // the composite quantity
        q ~ beta(1, 1); 
      
        // Induced priors on p11 and p10
        target += beta_lpdf(theta11 / q | 10, 2);
        target += beta_lpdf(theta10 / q | 2, 10);
      
        // Jacobian adjustment
        target += -2 * log(q);
      
        // Soft constraint: p11 > p10
        {
          real delta = p11 - p10;
          delta ~ normal(0.2, 0.5) T[0, ];
        }
      
        // Likelihood
        y ~ bernoulli(p_obs);
      }
      
      generated quantities {
        vector[N] log_lik;            // pointwise log-likelihood
        vector[N] Pr_z_eq_1;          // posterior Pr(z_i = 1 | y_i)
        array[N] int<lower=0, upper=1> ypred;

        for (n in 1:N) {
          // log-likelihood
          log_lik[n] = bernoulli_lpmf(y[n] | p_obs[n]);
      
          // posterior predictive
          ypred[n] = bernoulli_rng(p_obs[n]);
      
          // posterior Pr(z_i = 1 | y_i)
          {
            real log_prob_z1 = log(psi[n]);
            real log_prob_z0 = log1m(psi[n]);
      
            if (y[n] == 1) {
              log_prob_z1 += log(theta11);
              log_prob_z0 += log(theta10);
            } else {
              log_prob_z1 += log1m(theta11);
              log_prob_z0 += log1m(theta10);
            }
      
            Pr_z_eq_1[n] = inv_logit(log_prob_z1 - log_prob_z0);
          }
        }
      }
