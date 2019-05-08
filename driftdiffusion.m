function [i, restrictionPoint] = driftdiffusion(startI, generation)

thresholdX = 100;
thresholdY = 100;
checkpoint = 0;
diffuse_const=1;

xStart=0;

X(1,startI+1) = xStart;
Y(1,startI+1) = (0);
for k=1:startI
     X(k) = NaN;
     Y(k) = NaN;
end
i= startI + 1;
stop = 0;

percentG1 = 0.3;

if(generation == 0)
    %[myMin, age] = min(abs(rand(1)-invgcdf(0:0.1:50,1,0.5)));
    if(rand(1) > percentG1)
        checkpoint=1;
        restrictionPoint = i;
        Y=rand(1)*100;
    else
        X=rand(1)*100;
    end 
end


hold on;
while ~stop
    i=i+1;
    if checkpoint == 0
        X(i)= X(i-1) + diffuse_const + round(1*normrnd(0,1));
        Y(i)= Y(i-1) + round(normrnd(0,1));
    end
    if checkpoint == 1
        X(i)= X(i-1) + round(normrnd(0,1));
        Y(i)= Y(i-1) + diffuse_const + round(1*normrnd(0,1));
    end
    
    if checkpoint == 0 && X(i) > thresholdX
        checkpoint = 1;
        restrictionPoint = i;
        %Y(i)=0;
    end
    if checkpoint == 1 && Y(i) > thresholdY
       %checkpoint = 0;
       plot(i,Y(i),'X','MarkerSize', 10, 'linewidth', 2);
       stop = 1;
    end
end
plot(X,'color', [1-generation*0.2,0,0,0.2]);
plot(Y,'color', [0,1-generation*0.2,0,0.2]);
xEnd = X(end);
end
