function [g1_norm, g2_norm] = cell_pop_simp2(a1max,a2max,T,beta,init_type)
%B1 and B2 are the per capita rates of leaving a given stage.
%cells that leave the first stage enter the second stage.
%cells that leave the second stage, divide, and their daughters enter the
%first stage.

h=.01;

%gives num_steps age gridpoints for each age step when integrating along
%characteristics.
num_steps=3;
hfine=h/num_steps;

%get transition rates for each stage on a fine grid.
[beta1, beta2] = cell_pop_beta(a1max,a2max,T,hfine,beta);

%get initial population densities
[init_pop1, init_pop2] = cell_pop_initial_density(a1max,a2max,T,init_type);

%times for plotting;
t=0:h:T;
%ages for plotting;
ages1=0:h:a1max+T;
ages2=0:h:a2max+T;

[g1_norm, g2_norm, G1, G2] = cell_pop_PDE_solver(ages1,ages2,t,beta1,beta2,init_pop1, init_pop2,h,num_steps);


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

