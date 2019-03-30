function dbeta_plus_star = inverse_Gaussian_tp_nonneg_ODE(a,beta_plus_star,m,s)
%returns the rhs of the ODE that the inverse Gaussian tp satisfies

star=(-3./(2*a))+(1./(2*s^2*a.^2))-(m^2)/(2*s^2);
dstar=(3./(2*a.^2))-(1./(s^2*a.^3));

dbeta_plus_star = star.*(beta_plus_star-star)+(beta_plus_star-star).^2+dstar;

end