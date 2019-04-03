function [g1_norm, g2_norm] = cell_pop_simp2(a1max,a2max,T,beta)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.


Gaussian_density=@(a,mu,sigmasq)(1/(2*pi*sigmasq)^.5)*exp(-((a-mu).^2)/(2*sigmasq));

h=.01;
%gives num_steps age gridpoints for each age step when integrating along
%characteristics.
num_steps=3;
hfine=.01/num_steps;

%ages at which we approximate the population density
ages1=0:h:a1max+T;
ages2=0:h:a2max+T;

ages1_fine=0:hfine:a1max+T;
ages2_fine=0:hfine:a2max+T;

mu1=4;
mu2=2;
sigmasq1=.1;
sigmasq2=.1;

% Set the initial populations
init_pop1=Gaussian_density(ages1, mu1, sigmasq1)*10^(-2);
init_pop1(ages1>a1max)=0;
init_pop2=Gaussian_density(ages2,mu2, sigmasq2)*10^(-2);
init_pop2(ages2>a2max)=0;

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


if strcmp(beta,'gaussian')
% Transition probabilities have an age-dependency with Gaussian dist.

tp1=@(a)Gaussian_density(a,C1,1);
tp2=@(a)Gaussian_density(a,C2,1);

f1=@(a)B1*tp1(a);
f2=@(a)B2*tp2(a);

end

if strcmp(beta,'exponential')
tp1=@(a)ones(size(a));
tp2=@(a)ones(size(a));

f1=@(a)B1*tp1(a);
f2=@(a)B2*tp2(a);

end

if strcmp(beta,'inverse_gaussian')

f1=inverse_Gaussian_tp3(ages1_fine,.08,.02);
f2=inverse_Gaussian_tp3(ages2_fine,.08,.02);

f1_coarse=f1(mod(ages1_fine,3)==1);
f2_coarse=f2(mod(ages2_fine,3)==1);

end

%get the integral of the trasition rates over each age interval
int_f1=zeros(1,n1-1);
int_f2=zeros(1,n2-1);
for i1=1:n1-1
        f1_vector=f1((i1-1)*num_steps+1:i1*num_steps+1);
        int_f1(i1)=h_fine*trapz(f1_vector);
end

for i2=1:n2-1
        f2_vector=f2((i2-1)*num_steps+1:i2*num_steps+1);
        int_f2(i2)=h_fine*trapz(f2_vector);
end

%j for time
%i1 and i2 for ages
for j=1:m-1
    
    %get population density at next time step for next age by integrating
    %one step along the characterisitc.
    for i1=1:n1-1
        g1(i1+1,j+1)=g1(i1,j)*exp( -int_f1(i1) );
    end
    
    for i2=1:n2-1
        g2(i2+1,j+1)=g2(i2,j)*exp( -int_f2(i2) );
    end
    
    %get density of individuals of age zero, by integrating across all ages.
    %here the integral is approximated as a right-hand Riemann sum.
    G2(j+1)=h*sum(g2(2:n2,j+1));
    G1(j+1)=h*sum(g1(2:n1,j+1));
    
    % At the end of G2, two new cells enter G1 (cell division)
    %g1(1,j+1)=int_0^a1max(beta(a)g2(a,j)da)
    g1(1,j+1)=2*h*trapz(f1_coarse.*g2(:,j)) ;
    g2(1,j+1)=h*trapz(f2_coarse.*g1(:,j));
    
    g1_norm(:,j+1)=g1(:,j+1)/(G1(j+1)+g1(1,j+1)*h);
    g2_norm(:,j+1)=g2(:,j+1)/(G2(j+1)+g2(1,j+1)*h);
end

    hold off
    surf(t,ages1,g1_norm)
    hold on
    axis([0 T 0 T]);
    shading interp
    saveas(gcf,'G1')
    hold off
    
    surf(t,ages2,g2_norm)
    hold on
    axis([0 T 0 T]);
    shading interp
    saveas(gcf,'G2')
    
    hold off
    
    plot(t,G1./(G1+G2),'g','LineWidth',4);
    hold on 
    plot(t,G2./(G1+G2),'r','LineWidth',4);
    saveas(gcf,'fractions_G1_G2');

end

