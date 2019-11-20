function main()

mu1=.25;
sigma1=1;
mu2=.064;
sigma2=.031;

a1max=30;
a2max=30;
T=100;


[g, f, exp_uni_growth_rate_fit, exp_uni_mean1, exp_uni_var1, exp_uni_mean2, exp_uni_var2] = cell_pop_wrapper(a1max,a2max,T,'exponential','uniform',mu1,sigma1,mu2,sigma2);
[g, f, exp_gau_growth_rate_fit, exp_gau_mean1, exp_gau_var1, exp_gau_mean2, exp_gau_var2] = cell_pop_wrapper(a1max,a2max,T,'exponential','gaussian',mu1,sigma1,mu2,sigma2);
[g, f, exp_sta_growth_rate_fit, exp_sta_mean1, exp_sta_var1, exp_sta_mean2, exp_sta_var2] = cell_pop_wrapper(a1max,a2max,T,'exponential','stable_exp',mu1,sigma1,mu2,sigma2);
[g, f, invg_uni_growth_rate_fit, invg_uni_mean1, invg_uni_var1, invg_uni_mean2, invg_uni_var2] = cell_pop_wrapper(a1max,a2max,T,'inverse_gaussian','uniform',mu1,sigma1,mu2,sigma2);
[g, f, invg_gau_growth_rate_fit, invg_gau_mean1, invg_gau_var1, invg_gau_mean2, invg_gau_var2] = cell_pop_wrapper(a1max,a2max,T,'inverse_gaussian','gaussian',mu1,sigma1,mu2,sigma2);
[g, f, invg_sta_growth_rate_fit, invg_sta_mean1, invg_sta_var1, invg_sta_mean2, invg_sta_var2] = cell_pop_wrapper(a1max,a2max,T,'inverse_gaussian','stable_invg',mu1,sigma1,mu2,sigma2);


fprintf("%17s %19s %10s %10s %10s %10s %10s\n", "Model", "Initial Conditions", "Mean_g", "Variance_g", "Mean_f", "Variance_f", "Growth rate fit");
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Exponential", "Uniform", exp_uni_mean1, exp_uni_var1, exp_uni_mean2, exp_uni_var2, exp_uni_growth_rate_fit);
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Exponential", "Gaussian", exp_gau_mean1, exp_gau_var1, exp_gau_mean2, exp_gau_var2, exp_gau_growth_rate_fit);
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Exponential", "Stable", exp_sta_mean1, exp_sta_var1, exp_sta_mean2, exp_sta_var2, exp_sta_growth_rate_fit);
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Inverse Gaussian", "Uniform", invg_uni_mean1, invg_uni_var1, invg_uni_mean2, invg_uni_var2, invg_uni_growth_rate_fit);
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Inverse Gaussian", "Gaussian", invg_gau_mean1, invg_gau_var1, invg_gau_mean2, invg_gau_var2, invg_gau_growth_rate_fit);
fprintf("%17s %19s %1.4e %1.4e %1.4e %1.4e %1.4e\n", "Inverse Gaussian", "Stable", invg_sta_mean1, invg_sta_var1, invg_sta_mean2, invg_sta_var2, invg_sta_growth_rate_fit);
fprintf("\n");
convertfigs
end