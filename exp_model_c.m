% Compute c for the exponential model
function c = exp_model_c(beta_g, beta_f)
    c = (-(beta_g(1)+beta_f(1))+sqrt((beta_g(1)+beta_f(1))^2 + 4*beta_g(1)*beta_f(1)))/2;
    %c=0.045;
end