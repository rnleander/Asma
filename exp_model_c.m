% Compute c for the exponential model
function c = exp_model_c(beta_g, beta_f)
    c = (-(beta_g+beta_f)+sqrt((beta_g+beta_f)^2 + 4*beta_g*beta_f))/2;
end