function [y] = inverse_Gaussian_tp(a,m,s)
%computes the transition probability at the ages in a
a_fine=0:h:max(a);
y_fine=onestagepdf2(a_fine,m,s);
y1=onestagepdf1(a);
int_index=zeros(length(a));
int=zeros(length(a));
y=zeros(length(a));
for i=1:length(a)
int_index(i)=find(a_fine<=a(i),'last');
int(i)=h*sum(y_fine(0:int_index(i)));
y(i)=y1(i)/(1-int(i));
end


Simp_Rule(onestagepdf2,a,b,nstps)


end

