function dbeta = inverse_Gaussian_tp_ODE(a,beta,m,s)
%returns the rhs of the ODE that the inverse Gaussian tp satisfies

star=(-3./(2*a))+(1./(2*s^2*a.^2))-(m^2)/(2*s^2);

dbeta = star.*beta+beta.^2;

end