function mycell = experiment(startI, generation, genLimit,h,plotDDT,noGrowth)
    
mycell.generation = generation;
mycell.begin = startI;   
[i, restrictionPoint] = driftdiffusion(mycell.begin, mycell.generation,h,plotDDT);
mycell.imt = i - mycell.begin;
mycell.restrictionPoint = restrictionPoint;
mycell.end = i;

if generation < genLimit

    % run children
    leftdaughter = experiment(mycell.end, generation+1, genLimit,h,plotDDT,noGrowth);
    if(noGrowth==0)
        rightdaughter = experiment(mycell.end, generation+1, genLimit,h,plotDDT,noGrowth);
        mycell.progeny = {leftdaughter, rightdaughter};
    else
        mycell.progeny = {leftdaughter};
    end
else
    mycell.progeny = {};
end
end