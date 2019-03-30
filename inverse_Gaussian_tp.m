function [a,y,y1,denom,bb,lbt] = inverse_Gaussian_tp(a,m,s)
%computes the transition probability at the ages in a

%upper bound on beta
bb=(9/8)*s^2+m^2/(2*s^2);
%lower bound on t*
lbt=1/(3*s^2);
h=.0001;
a_fine=0:h:max(a);
y_fine=onestagepdf2(a_fine,m,s);
y1=onestagepdf2(a,m,s);
int_index=zeros(size(a));
int=zeros(size(a));
psum=zeros(size(a));
y=zeros(size(a));
for j=1:length(a)
int_index(j)=find(a_fine<=a(j),1,'last');
[psum(j),sum_failed] = sum2(y_fine(1:int_index(j)));
int(j)=psum(j)*h;
y(j)=y1(j)/(1-int(j));
end
denom=1-int;
end

