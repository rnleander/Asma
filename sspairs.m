function pairs = sspairs(ancestor)
    pairs = {};
    
    if size(ancestor.progeny,2) == 2
        parent = ancestor;
        pair.left = parent.progeny{1};
        pair.right = parent.progeny{2};
        
        leftpairs = sspairs(pair.left);
        rightpairs = sspairs(pair.right);
        
        pairs = cat(2,pair,leftpairs,rightpairs); 
    end
end