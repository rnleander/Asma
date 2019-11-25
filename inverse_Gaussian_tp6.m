function y = inverse_Gaussian_tp6(a,m,s, label,beta,init_type, perturb_str)
%approximates the transition probability at the ages in a using
%interpolation
%tic;
%clf;
set(figure, 'WindowStyle', 'docked')
hold on;
error_threshold = 0.1;
min_a = min(a);
max_a = max(a);
num_a = size(a,2);
h = 0.5*(min_a+max_a);
old_interpolant = rand(num_a,1)';
Legend=cell(1,1);
%fprintf("Starting interpolation with h=%f, min_a=%f, max_a=%f, num_a=%d, and error_threshold=%f\n", h, min_a, max_a, num_a, error_threshold);
idx=1;

% Uncomment this line to disable memoization between runs
clear sampleMemoize
num_knots={};
error={};
while true
    t = min_a:h:max_a;
    num_knots{idx} = size(t,2);
    if num_knots{idx} > num_a
        %fprintf("Interpolation failed, direct calculation instead of %d knots\n", num_knots);
        %tic;
        new_interpolant = inverse_Gaussian_tp3(a,m,s);
        %toc;
        break;
    end
    %fprintf("Evaluating %d knots with h=%f\n", num_knots, h);
    %tic;
    
    %k = inverse_Gaussian_tp3(t,m,s);
    
    % Uncomment this line to disable memoization entirely
    %clear sampleMemoize 
    k = sampleMemoize(t,m,s,h);
    
    %toc;
    new_interpolant = spline(t,k,a);
    %new_interpolant = pchip(t,k,a);
    error{idx} = norm(new_interpolant-old_interpolant);
    Legend{idx}=sprintf("n=%d error=%f", num_knots{idx}, error{idx});
    %fprintf("L2 distance between steps is %f\n\n", error);
    

    
    
    
    if error{idx} > error_threshold
       plot(a, new_interpolant, 'Color', [      0    0.4470    0.7410]);
       old_interpolant = new_interpolant;
       h=h/2;
       idx=idx+1;
    else
        plot(a, new_interpolant, 'Color', [0.8500    0.3250    0.0980]);
        %fprintf("Accepted interpolant with %d knots with %f error\n", num_knots, error);
        break
    end
    xlabel('Time (hrs)');
    ylabel('Maturation rate');
    title("Interpolation with n knots");
    title(sprintf("\\beta_%s", label));
    %legend(Legend);
    drawnow;
end

axh2 = get(gca, 'Children');
l=legend([axh2(1), axh2(2)], 'Accepted', 'Rejected');
if strcmp(label, "f") == 1
    p=l.Position;
    set(l, 'Position', [0.5, 0.4, p(3), p(4)]);
end

plot_filename = sprintf("figures/%s_%s_interpolation_beta_%s%s", beta, init_type, label, perturb_str);
saveas(gcf, plot_filename);

set(figure, 'WindowStyle', 'docked')
hold on;
knots_mat = cell2mat(num_knots);
error_mat = cell2mat(error);
plot(knots_mat(1:size(knots_mat,2)-1), error_mat(1:size(error_mat,2)-1),'*', 'Color', [      0    0.4470    0.7410])
plot(knots_mat(size(knots_mat,2)), error_mat(size(error_mat,2)),'*', 'Color', [0.8500    0.3250    0.0980])
plot_filename = sprintf("figures/%s_%s_interpolation_error_beta_%s%s", beta, init_type, label, perturb_str);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel("Number of knots");
ylabel(sprintf('\\boldmath${\\vert\\vert\\hat{\\beta}_{%s,i}-\\hat{\\beta}_{%s,i+1}\\vert\\vert}_2$',label, label),'Interpreter','latex')
title(sprintf("\\beta_%s", label));
axh2 = get(gca, 'Children');
% l=legend([axh2(1), axh2(2)], 'Accepted', 'Rejected');
% if strcmp(label, "f") == 1
%     p=l.Position;
%     set(l, 'Position', [0.3, 0.4, p(3), p(4)]);
% end
saveas(gcf, plot_filename);
y=new_interpolant;
%toc;
end
