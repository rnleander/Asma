function y = sampleMemoize(a, m, s, h)

%persistent last_a
persistent last_m
persistent last_s
persistent last_h
persistent last_y

y=zeros(1,size(a,2));
if isempty(last_h) || h ~= last_h/2
    %fprintf("Direct evaluation\n");
    y = inverse_Gaussian_tp3(a,m,s); 
elseif m == last_m && s == last_s
    %fprintf("Matching parameters\n");
    miss_idx=1;
    for a_idx=1:size(a,2)
        last_idx = (a_idx-1)/2+1;
        %fprintf("last_idx=%d a_idx=%d\n", last_idx, a_idx);           
        if floor(last_idx)==last_idx
            %fprintf("hit\n");
            %fprintf("last_a(last_idx)=%f a(a_idx)=%f\n", last_a(last_idx), a(a_idx));           
            y(a_idx) = last_y(last_idx);
        elseif a(a_idx)==0 % This is here because of a bug in tp3 for scalar 0
            fprintf("zero\n");
            y(a_idx) = 0;
        else
            %fprintf("miss\n");
            miss_a(miss_idx) = a(a_idx);
            miss_idx = miss_idx + 1;
            %y(a_idx) = inverse_Gaussian_tp3(a(a_idx), m, s);
        end 
    end
    miss_y = inverse_Gaussian_tp3(miss_a, m, s);
    miss_idx=1;
    for a_idx=1:size(a,2)
        last_idx = (a_idx-1)/2+1;
        %fprintf("last_idx=%d a_idx=%d\n", last_idx, a_idx);
        if floor(last_idx)~=last_idx
            %fprintf("miss two\n");
            y(a_idx) =  miss_y(miss_idx);
            miss_idx = miss_idx + 1;
        end
    end
end
%last_a = a;
last_m = m;
last_s = s;
last_h = h;
last_y = y;
end