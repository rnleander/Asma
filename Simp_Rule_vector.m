function [int_ab] = Simp_Rule_vector(f,h)
%f must have an even number of elements
assert(mod(length(f),2)==1,'f is not of even length\n');

ulsum=length(f)/2;

feven = f(2*(1:ulsum));
fodd = f(2*(1:ulsum-1)+1);
int_ab=(4*sum(feven)+2*sum(fodd)+f(1)+f(length(f)))*h/3;
end
