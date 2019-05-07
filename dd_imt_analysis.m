clear;
clf;
hold on;
frame(1:900)=250;
plot(frame);

M=[];
D=[];
L=[];
R=[];
for k=1:20
    ancestor = experiment(50,0,0,4);

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





