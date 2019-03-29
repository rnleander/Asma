function [y] = inverse_Gaussian_tp3(a,m,s)
%computes the transition probability at the ages in a

y1=onestagepdf2(a,m,s);

int = invgcdf(a,m,s);

y = y1./(1.-int);

end
