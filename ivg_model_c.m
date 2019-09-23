% find c for the inverse gaussian model
function c = ivg_model_c(lam_g, mu_g, lam_f, mu_f)

syms s

l_invg = @(lam, mu) sqrt(1/lam)*sqrt(lam)*exp(((lam)/(mu))-((sqrt(((lam)/(mu^2))+2*s))/(sqrt(1/lam))));

relation = (1/2) == l_invg(lam_g,mu_g)*l_invg(lam_f,mu_f);

c = solve(relation, s);

end