function [beta1, beta2] = cell_pop_beta(a1max,a2max,T,hfine,beta)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

%We will use parameters that correspond to the ML parameters for untreated 
%MCF10A cells. mu1=.25, sigma1=1, mu2=.064, sigma2=.031. Using the
%mean=1/mu and var=sigma^2/mu^3.

mu1=.25;
sigma1=1;
mu2=.064;
sigma2=.031;

mean1=1/mu1;
var1=sigma1^2/mu1^3;
mean2=1/mu2;
var2=sigma2^2/mu2^3;

ages1_fine=0:hfine:a1max+T;
ages2_fine=0:hfine:a2max+T;

%This option is currently nonfunctional
% if strcmp(beta,'gaussian')
% % Transition probabilities have an age-dependency with Gaussian dist.
%     mu=10;
%     sigmasq=1;
%     B1=10;
%     B2=100;
%     beta1=Gaussian_density(ages1_fine,mu,sigmasq);
%     beta2=Gaussian_density(ages2_fine,mu,sigmasq);
% 
%     beta1=B1*beta1;
%     beta2=B2*beta2;
% 
% end

if strcmp(beta,'exponential')
    lambda1=1/mean1;
    lambda2=1/mean2;
    beta1=ones(size(ages1_fine));
    beta2=ones(size(ages2_fine));

    beta1=lambda1*beta1;
    beta2=lambda2*beta2;

end

if strcmp(beta,'inverse_gaussian')

    beta1=inverse_Gaussian_tp6(ages1_fine,mu1,sigma1, "g");
    beta2=inverse_Gaussian_tp6(ages2_fine,mu2,sigma2, "f");


end

end
