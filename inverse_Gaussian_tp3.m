function y = inverse_Gaussian_tp3(a,m,s)
%computes the transition probability at the ages in a

a = vpa(a);

pdf = onestagepdf2(a,m,s);

cdf = invgcdf(a,m,s);

denom = 1-cdf;

y = pdf./(denom);

y = double(y);
end
