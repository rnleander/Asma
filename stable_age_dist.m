% find the stable age distribution
%
% currently this only handles inverse gaussian
% this needs to be extended to handle the exponential model as well
%
% TEST PROCEDURE:
% cell_pop_wrapper(30,30,100,'inverse_gaussian','gaussian')
% cell_pop_wrapper(30,30,100,'inverse_gaussian','stable_invg')
function g_vec_invgcdf = stable_age_dist(g_0, m_g, s_g, m_f, s_f, h, t_max)
%t_max = 50
%h = 0.1
%m_g = 0.1
%s_g = 0.08
%g_0 = 1

c = ivg_model_c(1/sqrt(m_g), 1/m_g, sqrt(1/m_f), 1/m_f)
s_mesh = 0:h:t_max;


figure
set(gcf, 'WindowStyle', 'docked')
hold on


% % This is a numerical method for finding the stable distribution
% g_numerical = @(idx, a, m_g, s_g, h, t_max) g_0*exp(-trapz(c + inverse_Gaussian_tp6(0:h:t_max,m_g,s_g)))
% g_vec_numerical=zeros(size(s_mesh,2),1);
% idx=1;
% for a=s_mesh
%     %idx = floor(a/h+1);
%     g_vec_numerical(idx) = g_numerical(idx, a, m_g, s_g, h, t_max);
%     idx = idx + 1;
% end
% trapz(g_vec_numerical)*h
% plot(s_mesh, g_vec_numerical)
% %plot(s_mesh, zero_patch_vector(g_vec_numerical))


% % This is an analytical method for finding the stable distribution
% g_dist = makedist('InverseGaussian','mu',1/m_g,'lambda',1/(s_g^2));
% g_analytic = @(a, c) g_0*(exp(-c*a)*(1-cdf(g_dist,a)))
% g_vec_analytic=zeros(size(s_mesh,2),1);
% idx=1;
% for a=s_mesh
%     %idx=floor(a/h+1)
%     %a = vpa(a);
%     %c = vpa(c);
%     g = g_analytic(a, c);
%     g_vec_analytic(idx) = g;
%     idx = idx + 1;
% end
% trapz(g_vec_analytic)*h
% plot(s_mesh, g_vec_analytic);
% %plot(s_mesh, zero_patch_vector(g_vec_analytic))


% Second analytic method
g_invgcdf = @(a, c, m, s) g_0*exp(-c*a)*(1-invgcdf(a, m, s))
g_vec_invgcdf=zeros(size(s_mesh,2),1);
idx=1;
for a=s_mesh
    g = g_invgcdf(a,c, m_g, s_g);
    if ~isequal(size(g),[1, 1])
        fprintf("OH NOES!\n");
    end   
    g_vec_invgcdf(idx) = g;
    idx = idx + 1;
end
trapz(g_vec_invgcdf)*h
plot(s_mesh, g_vec_invgcdf);
%plot(s_mesh, zero_patch_vector(g_vec_invgcdf))


hold off
end