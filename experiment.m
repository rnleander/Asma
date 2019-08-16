function mycell = experiment(startI, generation, genLimit,h,plotDDT,noGrowth,model,startingOffset)
    
mycell.generation = generation;
mycell.begin = startI;   
if(model==0)
[i, restrictionPoint] = driftdiffusion(mycell.begin, mycell.generation,h,plotDDT,startingOffset);
elseif(model==1)
[i, restrictionPoint] = exponential_process(mycell.begin, mycell.generation,h,plotDDT,startingOffset);
else
    disp("ERROR UNSUPPORTED MODEL")
end
mycell.imt = i - mycell.begin;
mycell.restrictionPoint = restrictionPoint;
mycell.end = i;

if generation < genLimit

    % run children
    leftdaughter = experiment(mycell.end, generation+1, genLimit,h,plotDDT,noGrowth,model,startingOffset);
    if(noGrowth==0)
        rightdaughter = experiment(mycell.end, generation+1, genLimit,h,plotDDT,noGrowth,model,startingOffset);
        mycell.progeny = {leftdaughter, rightdaughter};
    else
        mycell.progeny = {leftdaughter};
    end
else
    mycell.progeny = {};
end
end