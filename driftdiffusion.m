function [i, restrictionPoint] = driftdiffusion(startI,generation,h,plotDDT)

invg=0;

thresholdX = 1;
thresholdY = 1;
checkpoint = 0;

%mu1=.25;
%sigma1=1;
%mu2=.064;
%sigma2=.031;

diffuse_const1=1;
drift_const1=.25;
diffuse_const2=.031;
drift_const2=.064;

%diffuse_const1=1;
%drift_const1=1;
%diffuse_const2=1;
%drift_const2=1;

%parameters are such that time is measured in hours
%timestep=h 

%mu1=2;
%mu2=8;
%sigmasq1=8;
%sigmasq2=1.9;

avg_age1=2;
avg_age2=8;
stnd_dev1=8^.5;
stnd_dev2=1.9^.5;

if(invg==1)
    % Inverse Gaussian
    xStart=0;
else
    % Experimental Exponential distribution
    xStart=0;
end

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
    if(rand(1) > percentG1)
        age=normrnd(avg_age2,stnd_dev2);
        if age<0
            age=0;
        end
        checkpoint=1;
        restrictionPoint = i;
        Y(i)=normrnd(drift_const2*age,diffuse_const2*(age)^.5);
    else
        age=normrnd(avg_age1,stnd_dev1);
        if age<0
            age=0;
        end
        X(i)=normrnd(drift_const1*age,diffuse_const1*(age)^.5);
    end 
end

hold on;
if(invg==1)
    % Inverse Gaussian
    while ~stop
        i=i+1;
        if checkpoint == 0
            X(i)= X(i-1) + drift_const1*h + diffuse_const1*normrnd(0,h);
            Y(i)= NaN;
        end
        if checkpoint == 1
            X(i)= NaN;
            Y(i)= Y(i-1) + drift_const2*h + diffuse_const2*normrnd(0,h);
        end

        if checkpoint == 0 && X(i) > thresholdX
            checkpoint = 1;
            restrictionPoint = i;
            Y(i)=0;
        end
        if checkpoint == 1 && Y(i) > thresholdY
           %checkpoint = 0;
           %plot(i,Y(i),'X','MarkerSize', 10, 'linewidth', 2);
           stop = 1;
        end
    end

else
    
    % Experimental Exponential distribution
    while ~stop
        i=i+1;
        if checkpoint == 0
            X(i)= normrnd(0,1);
            Y(i)= NaN;
        end
        if checkpoint == 1
            X(i)= NaN;
            Y(i)= normrnd(0,1);
        end

        if checkpoint == 0 && X(i) > thresholdX
            checkpoint = 1;
            restrictionPoint = i;
            Y(i)=0;
        end
        if checkpoint == 1 && Y(i) > thresholdY
           %checkpoint = 0;
           %plot(i,Y(i),'X','MarkerSize', 10, 'linewidth', 2);
           stop = 1;
        end
    end
end


if(plotDDT~=0)
    plot(X,'color', [1,0,0,0.2]);
    plot(Y,'color', [0,1,0,0.2]);
    if(plotDDT==2)
        drawnow();
    end
end
xEnd = X(end);
end
