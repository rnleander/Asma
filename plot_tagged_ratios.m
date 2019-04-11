function plot_tagged_ratios()
filename = "2011-10-21 Well C06 tracking edited.csv"
A = readtable(filename);
h = 0.1;
idx = 1;

g1_total = 0;
g2_total = 0;
birth_g1 = zeros(120/h,120/h);

clf

hold on
for row=1:size(A,1)
    birthtime(row) = table2array(A(row,'Birthtimeh'));
    g1_time(row) = table2array(A(row,'G1_time'));
    g2_time(row) = table2array(A(row,'S_G2_M_time'));
    g1_start(row) = birthtime(row);
    g2_start(row) = birthtime(row)+g1_time(row);
    g1_end(row) = g2_start(row);
    g2_end(row) = g2_start(row)+g2_time(row);
    g1_total = g1_total + g1_time(row);
    g2_total = g2_total + g2_time(row);
    
    birthIdx = round(birthtime(row)/h)+1;
    g1Idx = round(g1_time(row)/h)+1;
    
    birth_g1(birthIdx,g1Idx) = birth_g1(birthIdx,g1Idx) + 1;
    
    X(row,1) = birthtime(row);
    X(row,2) = g1_time(row);
    X(row,3) = g2_time(row);
end
g1_average = g1_total/size(A,1);
g2_average = g2_total/size(A,1);

idx = kmeans(X,3);

hold on
view(40,35)
plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',12)
plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',12)
plot3(X(idx==3,1),X(idx==3,2),X(idx==3,3),'g.','MarkerSize',12)
title('3 Kmeans clusters in Birth time, G1 time and G2 time')
xlabel('Birth time')
ylabel('G1 time')
zlabel('G2 time')

return

scatter3(birthtime, g1_time,g2_time,'o')


return

heatmap(A,'Birthtimeh','S_G2_M_time')

return

birth_g1
hist3(birth_g1)
return

histogram(g2_time, 'BinWidth', 0.1)
title('Frequency distribution of G2 duration in FUCCI cells')
xlabel('G2 duration')
ylabel('Count')
return

for time=0:h:120
    t(idx)=time;
    g1_count = 0;
    g2_count = 0;
    for row=1:size(A,1)
        if(time>g1_start(row) && time+h<g1_end(row))
            g1_count = g1_count + 1;
        end
        if(time>g2_start(row) && time+h<g2_end(row))
            g2_count = g2_count + 1;
        end
    end
    g1(idx)=g1_count;
    g2(idx)=g2_count;
    idx=idx+1;
    hold off;
    plot(t,100*g2./(g1+g2));
    title('Percent of FUCCI cells in G2 phase');
    xlabel('Time (hrs)');
    ylabel('Percentage');
    drawnow;
end
end