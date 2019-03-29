function [psum,sum_failed] = sum2(y)
oldpsum=0;
psum=0;
k=1;
sum_failed=0;
while k<=length(y)
psum=oldpsum+y(k);
if psum==oldpsum && y(k)~=0
    sum_failed=1;
    psum2=0;
    for i=k:length(y)
        psum2=psum2+y(i);
    end
    k=length(y)+1;
    psum=psum+psum2;
else
    oldpsum=psum;
    k=k+1;
end
end
end

