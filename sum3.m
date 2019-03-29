function [psum,sum_failed] = sum3(y)
oldpsum=0;
psum=0;
sum_failed=0;
sorted_y=sort(y);
for k=1:length(y)
psum=oldpsum+sorted_y(k);
if psum==oldpsum && sorted_y(k)~=0
    sum_failed=1;
end
    oldpsum=psum;
end
end
