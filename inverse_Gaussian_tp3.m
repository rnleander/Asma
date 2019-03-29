function y = inverse_Gaussian_tp3(a,m,s)
%computes the transition probability at the ages in a

%y1=onestagepdf2(a,m,s);


%y = y1./(1.-invgcdf(a,m,s));


y = onestagepdf2(a,m,s)./(1-invgcdf(a,m,s));
end
