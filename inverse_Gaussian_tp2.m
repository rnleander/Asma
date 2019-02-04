function [a,tp,denom,sum_is_one] = inverse_Gaussian_tp2(a,m,s)
h=.0001;
a_fine=0:h:max(a);
y_fine=onestagepdf2(a_fine,m,s);
int_index=zeros(1,length(a));
sum_failed=zeros(1,length(a));
tp=zeros(1,length(a));
sum_is_one=zeros(1,length(a));
denom=zeros(1,length(a));

for i=1:length(a)
int_index(i)=find(a_fine<=a(i),1,'last');
if i==1
    s1=max(1-h*sum2(y_fine(1:int_index(i))),realmin);
    s2=max(realmin,1-h*sum2(y_fine(1:int_index(i)+1)));
    Fx_plus_h=-log(s2);
    Fx=-log(s1);
    tp(i)=(Fx_plus_h-Fx)/(h);
    denom(i)=s1;
elseif i==length(a)
        s1=max(realmin,1-h*sum2(y_fine(1:int_index(i))));
        s2=max(realmin,1-h*sum2(y_fine(1:int_index(i)-1)));
        Fx_minus_h=-log(s2);
        Fx=-log(s1);
        tp(i)=(Fx-Fx_minus_h)/(h);
        denom(i)=s1;
else
    [int1,sum_failed(i-1)]=sum2(y_fine(1:int_index(i)-1));
    s1=max(realmin,1-h*int1);
    if s1==realmin
        sum_is_one(i-1)=1;
    end
    s2=max(realmin,1-h*sum2(y_fine(1:int_index(i)+1)));
    Fx_plus_h=-log(s2);
    Fx_minus_h=-log(s1);
    tp(i)=(Fx_plus_h-Fx_minus_h)/(2*h);
    denom(i)=s2;
end
end
end