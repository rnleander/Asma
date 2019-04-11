function plot_tagged_ratios()
filename = "2011-10-21 Well C06 tracking edited.csv"
A = readtable(filename);
h = 0.1;
idx = 1;

for row=1:size(A,1)
    birthtime = table2array(A(row,'Birthtimeh'));
    g1_time = table2array(A(row,'G1_time'));
    g2_time = table2array(A(row,'S_G2_M_time'));
    g1_start(row) = birthtime;
    g2_start(row) = birthtime+g1_time;
    g1_end(row) = g2_start(row);
    g2_end(row) = g2_start(row)+g2_time;
end

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