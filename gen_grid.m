function [gp] = gen_grid(a,h,b)
n=1+(b-a)/h;
if mod(n,1)~=0
    fprintf('b-a is not a multiple of h.\n');
    return
end
gp=zeros(1,n);
gp(1)=a;
for i=1:n-1
    gp(i+1)=gp(i)+h;
end
end

