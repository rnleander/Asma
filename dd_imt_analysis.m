clear;
clf;
hold on;
frame(1:900)=250;
plot(frame);

M=[];
D=[];
L=[];
R=[];
ancestors = {};
for k=1:50
    ancestor = experiment(0,0,5);
    
    ancestors = cat(2,ancestors,{ancestor});

    pairs = mdpairs(ancestor);
    for i=1:size(pairs,2)
       m(i) = pairs(i).mother.imt;
       d(i) = pairs(i).daughter.imt;
    end
    M=cat(2,m,M);
    D=cat(2,d,D);

    pairs = sspairs(ancestor);
    for i=1:size(pairs,2)
       l(i) = pairs(i).left.imt;
       r(i) = pairs(i).right.imt;
    end
    L=cat(2,l,L);
    R=cat(2,r,R);
end

plot(M,D,'o','MarkerSize', 5, 'linewidth', 2, 'color', [0.5,0,0,0.5]);
mdcorr = corr(M',D','type','Spearman')

plot(L,R,'*','MarkerSize', 7, 'linewidth', 1, 'color', [0,0,0.5,0.5]);
sscorr = corr(L',R','type','Spearman')

allCells={};
for ancestorIdx=1:size(ancestors,2)
    descendents = flattenDescendents(ancestors{ancestorIdx});
    allCells = cat(2,allCells,descendents);
end

maxT = 0;
for cellIdx=1:size(allCells,2)
    thisCell = allCells{cellIdx};
    if(thisCell.end > maxT)
        maxT = thisCell.end;
    end
end

alive=[];
g1=[];
g2=[];

allCells;

numAllCells = size(allCells,2)
for i=1:maxT
    alive(i) = 0;
    g1(i) = 0;
    g2(i) = 0;
    for cellIdx=1:size(allCells,2)
        thisCell = allCells{cellIdx};
        
        if (i >= thisCell.begin && i < thisCell.end)
            alive(i) = alive(i) + 1;
        end
        if (i >= thisCell.begin && i < thisCell.restrictionPoint)
            g1(i) = g1(i) + 1;
        end
        if (i >= thisCell.restrictionPoint && i < thisCell.end)
            g2(i) = g2(i) + 1;
        end
    end
    
    
end

plot(200*g1./alive,'y')
plot(200*g2./alive,'b')

