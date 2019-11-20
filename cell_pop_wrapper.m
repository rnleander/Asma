function [g1_norm, g2_norm, growth_rate_fit, mean1, var1, mean2, var2] = cell_pop_wrapper(a1max,a2max,T,beta,init_type,mu1,sigma1,mu2,sigma2)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

close all;
set(gcf, 'WindowStyle', 'docked');

fprintf("model %s init_type %s\n", beta, init_type);

%G1G2ratio=1/4;
if strcmp(beta, 'inverse_gaussian')
    G1G2ratio=0.286475;
elseif strcmp(beta, 'exponential')
    G1G2ratio=0.434779;
else
    fprintf("ERROR: invalid distribution, please choose inverse_gaussian or exponential\n");
    return;
end

h=.01;

%gives num_steps age gridpoints for each age step when integrating along
%characteristics.
num_steps=3;
hfine=h/num_steps;

%get transition rates for each stage on a fine grid.
[beta1, beta2, mean1, var1, mean2, var2] = cell_pop_beta(a1max,a2max,T,hfine,beta,init_type,mu1,sigma1,mu2,sigma2);

%get initial population densities
[init_pop1, init_pop2] = cell_pop_initial_density(a1max,a2max,T,h,init_type,G1G2ratio,beta,mu1,sigma1,mu2,sigma2);

%times for plotting;
t=0:h:T;
%ages for plotting;
ages1=0:h:a1max+T;
ages2=0:h:a2max+T;

[g1_norm, g2_norm, G1, G2] = cell_pop_PDE_solver(ages1,ages2,t,beta1,beta2,init_pop1, init_pop2,h,num_steps);


%     hold off
%     surf(t,ages1,g1_norm)
%     hold on
%     axis([0 T 0 T]);
%     xlabel('t');
%     ylabel('age');
%     zlabel('g1_norm');
%     shading interp
%     saveas(gcf,'G1')
%     %save(G1.mat,gcf,'-v7.3')
%     hold off

    figure;
    set(gcf, 'WindowStyle', 'docked');
    %nexttile;
    hold off;
    [Xq,Yq] = meshgrid(0:T/1000:T);
    %Xq = 0:T/1000:T;
    %Yq = 0:(a1max+T)/1000:a1max+T;
    Vq = interp2(t,ages1,g1_norm,Xq,Yq);
    surf(Xq,Yq,Vq);
    hold on;
    axis([0 T 0 T/2]);
    xlabel('time (hrs)');
    ylabel('age (hrs)');
    zlabel('g age density');
    title('g age density through time');
    shading interp;
    view([-90, -90, 90]);
    %saveas(gcf,'G1');
    hold off;
    plot_filename = sprintf("figures/%s_%s_g_age_structure", beta, init_type);
    saveas(gcf, plot_filename);
      
%     surf(t,ages2,g2_norm)
%     hold on
%     axis([0 T 0 T]);
%     xlabel('t');
%     ylabel('age');
%     zlabel('g2_norm');
%     shading interp
%     saveas(gcf,'G2')
%     %save(G2.mat,gcf,'-v7.3')

    figure;
    set(gcf, 'WindowStyle', 'docked');
    %nexttile;
    [Xq,Yq] = meshgrid(0:T/1000:T);
    %Xq = 0:T/1000:T;
    %Yq = 0:(a2max+T)/1000:a2max+T;
    Vq = interp2(t,ages2,g2_norm,Xq,Yq);
    surf(Xq,Yq,Vq);
    hold on;
    axis([0 T 0 T/2]);
    xlabel('time (hrs)');
    ylabel('age (hrs)');
    zlabel('f age density');
    title('f age density through time');
    shading interp;
    view([-90, -90, 90]);
    %saveas(gcf,'G2');
    hold off;
    plot_filename = sprintf("figures/%s_%s_f_age_structure", beta, init_type);
    saveas(gcf, plot_filename);
    
    figure;
    set(gcf, 'WindowStyle', 'docked');
    %nexttile;
    hold off;
    plot(t,G1./(G1+G2));
    hold on ;
    plot(t,G2./(G1+G2));
    xlabel('Time (hrs)');
    ylabel('Norm');
    ylim([0, 1]);
    title('g and f fractions as a function of time');
    legend('g','f');
    plot_filename = sprintf("figures/%s_%s_gfratio", beta, init_type);
    saveas(gcf, plot_filename);
    %saveas(gcf,'fractions_G1_G2');

    lastg1 = G1(end);
    lastg2 = G2(end);
    ratio = lastg1/lastg2;
    
    fraction = lastg1/(lastg1+lastg2);
    
    fprintf("End of experiment g1/g2 ratio: %f\n", ratio);
    fprintf("End of experiment g1/(g1+g2) fraction: %f\n", fraction);
    
    total_population = G1+G2;
    
    expfit = fit(t(:), total_population(:), 'exp1');
    
    growth_rate_fit = expfit.b;
    fprintf("growth rate fit %f", growth_rate_fit);
    
    figure;
    set(gcf, 'WindowStyle', 'docked');
    plot(t, total_population);
    hold on;
    plot(t, feval(expfit, t), '--');
    if strcmp(beta, "exponential")
        if strcmp(init_type, "uniform")
            title("Exponential Uniform");
        end
        if strcmp(init_type, "gaussian")
            title("Exponential Gaussian");
        end
        if strcmp(init_type, "stable_exp")
            title("Exponential Stable");
        end
    end
    if strcmp(beta, "inverse_gaussian")
        if strcmp(init_type, "uniform")
            title("Inverse Gaussian Uniform");
        end
        if strcmp(init_type, "gaussian")
            title("Inverse Gaussian Gaussian");
        end
        if strcmp(init_type, "stable_exp")
            title("Inverse Gaussian Stable");
        end
    end
            
    %title('Population over time with exponential growth fit');
    xlabel('Time (hrs)');
    ylabel('Population');
    ylim([0 100]);
    legend('Experimental', 'Fit');
    plot_filename = sprintf("figures/%s_%s_growthrate", beta, init_type);
    saveas(gcf, plot_filename);
    
    fprintf("\n\n");
end

