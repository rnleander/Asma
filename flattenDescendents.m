function allCells = flattenDescendents(ancestor)
allCells = {};
numChildren = size(ancestor.progeny,2);
allCells = cat(2, allCells, {ancestor} );
for childIdx=1:numChildren
    thisChild = ancestor.progeny{childIdx};
    descendents = flattenDescendents( thisChild );
    allCells = cat(2, allCells, descendents ); 
end
end