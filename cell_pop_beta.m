function [beta1, beta2, mean1, var1, mean2, var2] = cell_pop_beta(a1max,a2max,T,hfine,beta,init_type,mu1,sigma1,mu2,sigma2, perturb_str)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

%We will use parameters that correspond to the ML parameters for untreated 
%MCF10A cells. mu1=.25, sigma1=1, mu2=.064, sigma2=.031. Using the
%mean=1/mu and var=sigma^2/mu^3.

% mu1=.25;
% sigma1=1;
% mu2=.064;
% sigma2=.031;

mean1=1/mu1;
var1=sigma1^2/mu1^3;
mean2=1/mu2;
var2=sigma2^2/mu2^3;

ages1_fine=0:hfine:a1max+T;
ages2_fine=0:hfine:a2max+T;
ages_imt=0:hfine:a1max+T+a2max+T;

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
    
    var1=1/(lambda1^2);
    var2=1/(lambda2^2);
    
    imt1=lambda1*exp(-lambda1*ages1_fine);
    imt2=lambda2*exp(-lambda2*ages2_fine);
    
    beta1=ones(size(ages1_fine));
    beta2=ones(size(ages2_fine));

    beta1=lambda1*beta1;
    beta2=lambda2*beta2;
    
    label="Exponential";
    params1 = sprintf("\\lambda_g=%f", lambda1);
    params2 = sprintf("\\lambda_f=%f", lambda2);
end

if strcmp(beta,'inverse_gaussian')
    
    
    imt1=onestagepdf2(ages1_fine,mu1,sigma1);
    imt2=onestagepdf2(ages2_fine,mu2,sigma2);
    
    beta1=inverse_Gaussian_tp6(ages1_fine,mu1,sigma1, "g",beta,init_type, perturb_str);
    beta2=inverse_Gaussian_tp6(ages2_fine,mu2,sigma2, "f",beta,init_type, perturb_str);

    label="Inverse Gaussian";
    params1 = sprintf("\\mu_g=%f \\sigma_g=%f", mu1, sigma1);
    params2 = sprintf("\\mu_f=%f \\sigma_f=%f", mu2, sigma2);
end


imt1=imt1/(trapz(imt1)*hfine);
imt2=imt2/(trapz(imt2)*hfine);
imt=conv(imt1,imt2);
imt=imt/(trapz(imt)*hfine);

reduced_idx=ceil(50/hfine);
imt1=imt1(1:1:reduced_idx);
imt2=imt2(1:1:reduced_idx);
imt=imt(1:1:reduced_idx);
reduced_ages=ages_imt(1:1:reduced_idx);

figure;
set(gcf, 'WindowStyle', 'docked');
plot(reduced_ages, imt1);
hold on;
plot(reduced_ages, imt2);
plot(reduced_ages, imt);
%title(sprintf("Stage maturation time and total intermitotic time PDFs (%s)", label));
title(label);
xlabel('Time (hrs)');
ylabel('Probability Density');
%legend(sprintf("I_g(a); %s", params1), sprintf("I_f(a); %s", params2), "IMT(a)=I_g(a)*I_f(a)");
legend("I_g(a)", "I_f(a)", "(I_g \ast I_f)(a)");
drawnow;

plot_filename = sprintf("figures/%s_%s_imt%s", beta, init_type, perturb_str);
saveas(gcf, plot_filename);

end
