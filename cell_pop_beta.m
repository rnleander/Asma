function [beta1, beta2] = cell_pop_beta(a1max,a2max,T,hfine,beta)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

ages1_fine=0:hfine:a1max+T;
ages2_fine=0:hfine:a2max+T;

if strcmp(beta,'gaussian')
% Transition probabilities have an age-dependency with Gaussian dist.
    mu=10;
    sigmasq=1;
    B1=10;
    B2=100;
    beta1=Gaussian_density(ages1_fine,mu,sigmasq);
    beta2=Gaussian_density(ages2_fine,mu,sigmasq);

    beta1=B1*beta1;
    beta2=B2*beta2;

end

if strcmp(beta,'exponential')
    B1=10;
    B2=5;
    beta1=ones(size(ages1_fine));
    beta2=ones(ages2_fine);

    beta1=B1*beta1;
    beta2=B2*beta2;

end

if strcmp(beta,'inverse_gaussian')

    beta1=inverse_Gaussian_tp3(ages1_fine,.08,.02);
    beta2=inverse_Gaussian_tp3(ages2_fine,.08,.02);


end

end
