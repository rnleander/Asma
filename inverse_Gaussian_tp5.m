function [a,beta] = inverse_Gaussian_tp5(a,m,s)
hold off
%computes the transition probability at the ages in a
star=@(a)(-3./(2*a))+(1./(2*s^2*a.^2))-(m^2)/(2*s^2);
beta = onestagepdf2(a,m,s)./(1-invgcdf(a,m,s));
plot(a,beta,'b')
hold on
jac=@(t,y)J(t,y,m,s);
%opts = odeset('Stats','on','Jacobian',jac,'RelTol',1e-3,'AbsTol',1e-3);
%opts = odeset('Stats','on','Jacobian',jac,'NonNegative',1);
opts = odeset('Stats','on','Jacobian',jac,'NonNegative',1,'RelTol',1e-10,'AbsTol',1e-10,'MaxStep',1e-2);
[t,beta2_plus_star]=ode45(@(t,y)inverse_Gaussian_tp_nonneg_ODE(t,y,m,s),a(1600:length(a)),beta(1600)+star(a(1600)),opts);

plot(a(1600:length(a)),-star(a(1600:length(a))),'r');
plot(t,beta2_plus_star-star(t),'g');
beta(1600:length(a))=beta2_plus_star-star(a(1600:length(a)))';
hold off
end