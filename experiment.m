function mycell = experiment(startI, generation, genLimit)
    
mycell.generation = generation;
mycell.begin = startI;   
[i, restrictionPoint] = driftdiffusion(mycell.begin, mycell.generation);
mycell.imt = i - mycell.begin;
mycell.restrictionPoint = restrictionPoint;
mycell.end = i;

if generation < genLimit

    % run children
    leftdaughter = experiment(mycell.end, generation+1, genLimit);
    rightdaughter = experiment(mycell.end, generation+1, genLimit);
    
    mycell.progeny = {leftdaughter, rightdaughter};
else
    mycell.progeny = {};
end
end