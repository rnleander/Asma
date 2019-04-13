function [init_pop1, init_pop2] = cell_pop_initial_density(a1max,a2max,T,h,init_type,G1G2ratio)
%the parameters for the initial distribution are chosen so as to be
%reasonable given the distribution of times spent in the predicted first and second
%stages of the MCF10A data. 

%ages at which we approximate the population density
ages1=0:h:a1max+T;
ages2=0:h:a2max+T;

if strcmp(init_type,'gaussian')
Gaussian_density=@(a,mu,sigmasq)(1/(2*pi*sigmasq)^.5)*exp(-((a-mu).^2)/(2*sigmasq));

mu1=2;
mu2=8;
sigmasq1=8;
sigmasq2=1.9;

% Set the initial populations
init_pop1=Gaussian_density(ages1, mu1, sigmasq1)*10^(-2);
init_pop1(ages1>a1max)=0;
G1tot = cdf('Normal',a1max,mu1,sigmasq1^.5); 

init_pop2=Gaussian_density(ages2,mu2, sigmasq2)*10^(-2);
init_pop2(ages2>a2max)=0;
G2tot = cdf('Normal',a2max,mu2,sigmasq2^.5); 

init_pop1=init_pop1*G1G2ratio*G2tot/G1tot;

end

if strcmp(init_type,'uniform')
    
    init_pop1=zeros(size(ages1));
    init_pop2=zeros(size(ages2));

    pen_factor1=.5;
    radius_factor1=.5;
    pen_factor2=.5;
    radius_factor2=.5;

% Set the initial populations
    init_pop1_ub=pen_factor1*a1max+radius_factor1*min(a1max-pen_factor1*a1max,pen_factor1*a1max);
    init_pop1_lb=pen_factor1*a1max-radius_factor1*min(a1max-pen_factor1*a1max,pen_factor1*a1max);
    range1=init_pop1_ub-init_pop1_lb;

    init_pop1(init_pop1_lb <= ages1.*init_pop1_ub >= ages1)=1/range1;
    
    init_pop2_ub=pen_factor2*a2max+radius_factor2*min(a2max-pen_factor2*a2max,pen_factor2*a2max);
    init_pop2_lb=pen_factor2*a2max-radius_factor2*min(a2max-pen_factor2*a2max,pen_factor2*a2max);
    range2=init_pop2_ub-init_pop2_lb;

    init_pop2(init_pop2_lb <= ages2.*init_pop2_ub >= ages1)=1/range2;
    

end

end