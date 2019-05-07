function pairs = mdpairs(ancestor)
    pairs = {};
    
    if size(ancestor.progeny,2) == 2
        parent = ancestor;
        leftchild = parent.progeny{1};
        rightchild = parent.progeny{2};
      
        leftpair.mother= ancestor;
        leftpair.daughter=leftchild;
        
        rightpair.mother=ancestor;
        rightpair.daughter=rightchild; 
        
        leftpairs = mdpairs(leftchild);
        rightpairs = mdpairs(rightchild);
        
        pairs = cat(2,leftpair,rightpair,leftpairs,rightpairs);
    end
end