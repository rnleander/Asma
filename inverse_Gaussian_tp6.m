function y = inverse_Gaussian_tp6(a,m,s)
%approximates the transition probability at the ages in a using
%interpolation
tic;
clf;
hold on;
error_threshold = 0.001;
min_a = min(a);
max_a = max(a);
num_a = size(a,2);
h = 0.5*(min_a+max_a);
old_interpolant = rand(num_a,1)';
Legend=cell(1,1);
fprintf("Starting interpolation with h=%f, min_a=%f, max_a=%f, num_a=%d, and error_threshold=%f\n", h, min_a, max_a, num_a, error_threshold);
idx=1;

% Uncomment this line to disable memoization between runs
clear sampleMemoize

while true
    t = min_a:h:max_a;
    num_knots = size(t,2);
    if num_knots > num_a
        fprintf("Interpolation failed, direct calculation instead of %d knots\n", num_knots);
        tic;
        new_interpolant = inverse_Gaussian_tp3(a,m,s);
        toc;
        break;
    end
    fprintf("Evaluating %d knots with h=%f\n", num_knots, h);
    tic;
    
    %k = inverse_Gaussian_tp3(t,m,s);
    
    % Uncomment this line to disable memoization entirely
    %clear sampleMemoize 
    k = sampleMemoize(t,m,s,h);
    
    toc;
    new_interpolant = spline(t,k,a);
    %new_interpolant = pchip(t,k,a);
    error = norm(new_interpolant-old_interpolant);
    Legend{idx}=sprintf("n=%d error=%f", num_knots, error);
    fprintf("L2 distance between steps is %f\n\n", error);
    plot(a, new_interpolant);
    title("Interpolation with n knots");
    legend(Legend);
    drawnow;
    if error > error_threshold
       old_interpolant = new_interpolant;
       h=h/2;
       idx=idx+1;
    else
        fprintf("Accepted interpolant with %d knots with %f error\n", num_knots, error);
        break
    end 
end
y=new_interpolant;
toc;
end
