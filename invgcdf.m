function y = invgcdf(t,m,s)

mu = 1/m;

lam = 1/(s^2);

arg1 = sqrt(lam./t).*(t./mu - 1);

arg2 = -1.*sqrt(lam./t).*(t./mu + 1);

norm1 = cdf('Normal',arg1,0,1);

norm2 = cdf('Normal',arg2,0,1);

y = norm1+exp(2*lam/mu)*norm2;

end