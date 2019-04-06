function [g1_norm, g2_norm, G1, G2] = cell_pop_PDE_solver(ages1,ages2,t,beta1,beta2,init_pop1,init_pop2,h,num_steps)

%gives num_steps age gridpoints for each age step when integrating along
%characteristics.
hfine=h/num_steps;

%get transition rates on a coarser subgrid over which the PDE is solved.
beta1_coarse=beta1(mod(1:1:length(beta1),num_steps)==1);
beta2_coarse=beta2(mod(1:1:length(beta2),num_steps)==1);

n1=length(ages1);
n2=length(ages2);

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

%get the integral of the trasition rates over each age interval
int_beta1=zeros(1,n1-1);
int_beta2=zeros(1,n2-1);
for i1=1:n1-1
        beta1_vector=beta1((i1-1)*num_steps+1:i1*num_steps+1);
        int_beta1(i1)=hfine*trapz(beta1_vector);
end

for i2=1:n2-1
        beta2_vector=beta2((i2-1)*num_steps+1:i2*num_steps+1);
        int_beta2(i2)=hfine*trapz(beta2_vector);
end

%j for time
%i1 and i2 for ages
for j=1:m-1
    
    %get population density at next time step for next age by integrating
    %one step along the characterisitc.
    for i1=1:n1-1
        g1(i1+1,j+1)=g1(i1,j)*exp( -int_beta1(i1) );
    end
    
    for i2=1:n2-1
        g2(i2+1,j+1)=g2(i2,j)*exp( -int_beta2(i2) );
    end
    
    %get density of individuals of age zero, by integrating across all ages.
    %here the integral is approximated as a right-hand Riemann sum.
    G2(j+1)=h*sum(g2(2:n2,j+1));
    G1(j+1)=h*sum(g1(2:n1,j+1));
    
    % At the end of G2, two new cells enter G1 (cell division)
    %g1(1,j+1)=int_0^a1max(beta(a)g2(a,j)da)
    g1(1,j+1)=2*h*trapz(beta1_coarse'.*g2(:,j)) ;
    g2(1,j+1)=h*trapz(beta2_coarse'.*g1(:,j));
    
    g1_norm(:,j+1)=g1(:,j+1)/(G1(j+1)+g1(1,j+1)*h);
    g2_norm(:,j+1)=g2(:,j+1)/(G2(j+1)+g2(1,j+1)*h);
end