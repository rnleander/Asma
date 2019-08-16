function [i, restrictionPoint] = exponential_process(startI,generation,h,plotDDT,startingOffset)

checkpoint = 0;

avg_age1=2;
avg_age2=8;
stnd_dev1=8^.5;
stnd_dev2=1.9^.5;

theta = 500;

percentG1 = 0.3;

i= startI + 1;
age = 0;

if(generation == 0)
    if(rand(1) > percentG1)
        age=normrnd(avg_age2,stnd_dev2);

        checkpoint=1;
        restrictionPoint = i;
    else
        age=normrnd(avg_age1,stnd_dev1);

    end 
    if age<-startingOffset
        age=-startingOffset;
    end
end

if(checkpoint==0)
    restrictionPoint = startI + exprnd(theta) - age;
end
i = restrictionPoint + exprnd(theta) - age;


restrictionPoint = ceil(restrictionPoint);
i = ceil(i);

% if(restrictionPoint<1)
%     restrictionPoint=1;
% end
% if(i<=restrictionPoint)
%     i=restrictionPoint+1;
% end

end