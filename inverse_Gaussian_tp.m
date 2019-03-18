function [a,y,y1,denom] = inverse_Gaussian_tp(a,m,s)
%computes the transition probability at the ages in a
h=.001;
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
y(j)=1/((1/y1(j))-(int(j)/y1(j)));
end
denom=1-int;
end

