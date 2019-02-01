function [g1_norm, g2_norm] = cell_pop_simp2(B1,B2,C1,C2,T)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

eps=.01;
a1max=C1*10;
a2max=C2*10;


Gaussian_density=@(a,mu,sigmasq)(1/(2*pi*sigmasq)^.5)*exp(-((a-mu).^2)/(2*sigmasq));

h=.01;

%ages at which we approximate the population density
ages1=0:h:a1max;
ages2=0:h:a2max;

% Set the initial populations
init_pop1=Gaussian_density(ages1, 1, .1)*10^(-2);
init_pop2=Gaussian_density(ages2, 4, .1)*10^(-2);

n1=length(ages1);
n2=length(ages2);

t=0:h:T;

m=length(t);

%population density in each class
g1=zeros(n1,m);
g2=zeros(n2,m);

%normalized population density in each class (if we integrate g1_norm over
%all ages we get 1.)
g1_norm=zeros(n1,m);
g2_norm=zeros(n2,m);

%total number in nonzero age class (the integral of g1 and g2 across all nonzero ages.)
G1=zeros(1,m);
G2=zeros(1,m);

g1(:,1)=init_pop1;
g2(:,1)=init_pop2;

G2(1)=h*sum(g2(1:n2,1));
G1(1)=h*sum(g1(1:n1,1));

g1_norm(:,1)=g1(:,1)/G1(1);
g2_norm(:,1)=g2(:,1)/G2(1);

% Transition probabilities have an age-dependency with Gaussian dist.
%tp1=@(a)Gaussian_density(a,C1,1);
%tp2=@(a)Gaussian_density(a,C2,1);

tp1=@(a)ones(size(a));
tp2=@(a)ones(size(a));

f1=@(a)B1*tp1(a);
f2=@(a)B2*tp2(a);

%j for time
%i1 and i2 for ages
for j=1:m-1
    
    %get population density at next time step for next age by integrating
    %one step along the characterisitc.
    for i1=1:n1-1
        g1(i1+1,j+1)=g1(i1,j)*exp( -Simp_Rule(f1, 0, h, 10) );
    end
    
    for i2=1:n2-1
        g2(i2+1,j+1)=g2(i2,j)*exp( -Simp_Rule(f2, 0, h, 10) );
    end
    
    %get density of individuals of age zero, by integrating across all ages.
    %here the integral is approximated as a right-hand Riemann sum.
    G2(j+1)=h*sum(g2(2:n2,j+1));
    G1(j+1)=h*sum(g1(2:n1,j+1));
    
    % At the end of G2, two new cells enter G1 (cell division)
    g1(1,j+1)=2*B2*G2(j+1);
    g2(1,j+1)=B1*G1(j+1);
    
    g1_norm(:,j+1)=g1(:,j+1)/(G1(j+1)+g1(1,j+1)*h);
    g2_norm(:,j+1)=g2(:,j+1)/(G2(j+1)+g2(1,j+1)*h);
end

    hold off
    surf(t,ages1,g1_norm)
    hold on
    axis([0 T 0 T])
    shading interp
    saveas(gcf,'G1')
    hold off
    
    surf(t,ages2,g2_norm)
    hold on
    axis([0 T 0 T])
    shading interp
    saveas(gcf,'G2')

end

