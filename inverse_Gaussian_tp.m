function [a,y,y1,denom] = inverse_Gaussian_tp(a,m,s)
%computes the transition probability at the ages in a
h=.0001;
a_fine=0:h:max(a);
y_fine=onestagepdf2(a_fine,m,s);
y1=onestagepdf2(a,m,s);
int_index=zeros(size(a));
int=zeros(size(a));
y=zeros(size(a));
for j=1:length(a)
int_index(j)=find(a_fine<=a(j),1,'last');
oldpsum=0;
k=1;
while k<=int_index(j)
psum=oldpsum+y_fine(k);
if psum==oldpsum && y_fine(k)~=0
    psum2=0;
    for i=k:int_index(j)
        psum2=psum2+y_fine(i);
    end
    k=int_index(j)+1;
    psum=psum+psum2;
else
    oldpsum=psum;
    k=k+1;
end

end
int(j)=psum*h;
y(j)=y1(j)/(max(realmin,1-int(j)));
end
denom=1-int;
end

