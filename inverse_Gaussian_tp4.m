function [a,beta] = inverse_Gaussian_tp4(a,m,s)
%computes the transition probability at the ages in a
star=(-3./(2*a))+(1./(2*s^2*a.^2))-(m^2)/(2*s^2);
p=[-m^2/(2*s^2) -3/2 1/(2*s^2)];
r=roots(p);
r1=r(r>0);
r2=(-p(2)-(p(2)^2-4*p(1)*p(3))^.5)/(2*p(1));
beta=zeros(length(a),1);
beta(1:1600) = onestagepdf2(a(1:1600),m,s)./(1-invgcdf(a(1:1600),m,s));
jac=@(t,y)J(t,y,m,s);
%opts = odeset('Stats','on','Jacobian',jac,'RelTol',1e-3,'AbsTol',1e-3);
opts = odeset('Stats','on','Jacobian',jac);
%opts = odeset('Stats','on');
[t,beta2]=ode23s(@(t,y)inverse_Gaussian_tp_ODE(t,y,m,s),a(1600:length(a)),beta(1600),opts);
%[t,beta2]=ode45(@(t,y)inverse_Gaussian_tp_ODE(t,y,m,s),a(2200:length(a)),beta(2200),opts);
beta(1600:length(beta))=beta2;
end