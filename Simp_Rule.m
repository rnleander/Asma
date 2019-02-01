function [int_ab] = Simp_Rule(fun,a,b,nstps)
%nsteps=number of steps for integration.  It must be even
gridpts=linspace(a,b,nstps+1);
h=gridpts(2)-gridpts(1);
ulsum=nstps/2;

f=fun(gridpts);

feven = f(2*(1:ulsum));
fodd = f(2*(1:ulsum-1)+1);
int_ab=(4*sum(feven)+2*sum(fodd)+f(1)+f(nstps+1))*h/3;
end

