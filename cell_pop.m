function [] = cell_pop(B1,B2,T)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.


%This code is full of errors, cell_pop_simp contains the needed
%corrections.
eps=.01;
a1max=log(1/eps)/B1;
a2max=log(1/eps)/B2;

mu1=1/B1;
sigma1sq=.1;

mu2=1/B2;
sigma2sq=.1;

initial_pop_density=@(a,mu,sigmasq)(1/(2*pi*sigmasq)^.5)*exp(-((a-mu1).^2)/(2*sigma1sq));

h=.1;

%ages at which we approximate the population density
ages1=0:h:a1max;
ages2=0:h:a2max;

init_pop1=initial_pop_density(ages1,mu1,sigma1sq)*10^(-2);
init_pop2=initial_pop_density(ages2,mu2,sigma2sq)*10^(-2);

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
%total number in each class (the integral of g1 and g2 across all ages.)
G1=zeros(1,m);
G2=zeros(1,m);

g1(:,1)=init_pop1;
g2(:,1)=init_pop2;

G1(1)=h*sum(g2(2:n2,1));
G2(1)=h*sum(g1(2:n1,1));


 g1_norm(:,1)=g1(:,1)/G1(1);
 g2_norm(:,1)=g2(:,1)/G2(1);

%j for time
%i1 and i2 for ages
f1=@(t,y1)-B1*y1;
f2=@(t,y2)-B2*y2;
for j=1:m-1
    %get population density at next time step for next age by integrating
    %one step along the characterisitc.
    for i1=1:n1-1
        [T,Y1]=ode45(@(t,y1)f1(t,y1),[0,h],g1(i1,j));
        g1(i1+1,j+1)=Y1(2);
    end
    for i2=1:n2-1
        [T,Y2]=ode45(@(t,y2)f2(t,y2),[0,h],g2(i2,j));
        g2(i2+1,j+1)=Y2(2);
    end
    %get density of individuals of age zero, by integrating across all ages.
    %here the integral is approximated as a right-hand Riemann sum.
    G2(j+1)=h*sum(g2(2:n2,j+1));
    G1(j+1)=h*sum(g1(2:n1,j+1));
    g1(1,j+1)=2*B2*G2(j+1);
    g2(1,j+1)=B1*G1(j+1);
    
    g1_norm(:,j+1)=g1(:,j+1)/G1(j+1);
    g2_norm(:,j+1)=g2(:,j+1)/G2(j+1);
end


    
    surf(t,ages1,g1_norm)
    saveas(gcf,'G1')
    
    surf(t,ages2,g2_norm)
    saveas(gcf,'G2')
    
    
        









end

