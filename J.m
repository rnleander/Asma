function jac = J(a,beta,m,s)
%returns the jacobian of the rhs of the ODE that the inverse Gaussian tp satisfies
star=-(3/2)*a.^(-1)+(1/(2*s^2))*a.^(-2)-m/(2*s^2);
jac = star+2*beta;
end