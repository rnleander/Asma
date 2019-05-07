function mycell = experiment(inheritedX, startI, generation, genLimit)
    
mycell.generation = generation;
mycell.begin = startI;   
[i, X] = driftdiffusion(inheritedX, mycell.begin, mycell.generation);
mycell.xEnd = X;
mycell.imt = i - mycell.begin;
mycell.end = i;

if generation < genLimit
    % partition error
    leftInheritX = mycell.xEnd*0.6;
    rightInheritX = mycell.xEnd*0.4;
    
    % run children
    leftdaughter = experiment(leftInheritX, mycell.end, generation+1, genLimit);
    rightdaughter = experiment(rightInheritX, mycell.end, generation+1, genLimit);
    
    mycell.progeny = {leftdaughter, rightdaughter};
else
    mycell.progeny = {};
end
end