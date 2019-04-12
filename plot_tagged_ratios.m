function plot_tagged_ratios()
filename = "2011-10-21 Well C06 tracking edited.csv"
A = readtable(filename);
h = 0.1;



birth_g1 = zeros(120/h,120/h);

clf

maxtime = 0;
for row=1:size(A,1)
    death_time = table2array(A(row,'Birthtimeh')) + table2array(A(row,'G1_time')) + table2array(A(row,'S_G2_M_time'));
    if(death_time > maxtime)
        maxtime = death_time;
    end
end

hold on
censored_row=1;
g1_total = 0;
g2_total = 0;
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
    
    life_time(row) = g1_time(row);
    remaining_time(row) = maxtime-birthtime(row);
    
    lifetime_deficit(row) = remaining_time(row)-life_time(row);
    
    X(row,4) = lifetime_deficit(row);
    
    X(row,5) = g1_time(row)/g2_time(row);
    
    Y(row,1) = X(row,4);
    Y(row,3) = X(row,5);
    
    if(lifetime_deficit(row) > 0.3)
        X_censored(censored_row,1) = X(row,1);
        X_censored(censored_row,2) = X(row,2);
        X_censored(censored_row,3) = X(row,3);
        g1_start_censored(censored_row) = birthtime(row);
        g2_start_censored(censored_row) = birthtime(row)+g1_time(row);
        g1_end_censored(censored_row) = g2_start(row);
        g2_end_censored(censored_row) = g2_start(row)+g2_time(row);
    
        censored_row = censored_row + 1;
    end
    
end
g1_average = g1_total/size(A,1);
g2_average = g2_total/size(A,1);



hold on
view(40,35)

% idx = kmeans(X,3);
% plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r.','MarkerSize',12)
% plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.','MarkerSize',12)
% plot3(X(idx==3,1),X(idx==3,2),X(idx==3,3),'g.','MarkerSize',12)
% title('3 Kmeans clusters in Birth time, G1 time and G2 time')
% xlabel('Birth time')
% ylabel('G1 time')
% zlabel('G2 time')

% idx2 = kmeans(Y,3);
% plot(X(idx2==1,4),X(idx2==1,5),'r.','MarkerSize',12)
% plot(X(idx2==2,4),X(idx2==2,5),'b.','MarkerSize',12)
% title('3 Kmeans clusters in Birth time, G1 time, G2 time, Life time deficit and G1/G2 ratio')
% xlabel('Life time deficit')
% ylabel('G1/G2 ratio')

censored_idx = kmeans(X_censored,2);
% plot3(X_censored(censored_idx==1,1),X_censored(censored_idx==1,2),X_censored(censored_idx==1,3),'r.','MarkerSize',12)
% plot3(X_censored(censored_idx==2,1),X_censored(censored_idx==2,2),X_censored(censored_idx==2,3),'b.','MarkerSize',12)
% title('3 Kmeans clusters in Birth time, G1 time and G2 time after censoring')
% xlabel('Birth time')
% ylabel('G1 time')
% zlabel('G2 time')
% return

% scatter3(birthtime, g1_time,g2_time,'o')
% return

% heatmap(A,'Birthtimeh','S_G2_M_time')
% return

% birth_g1
% hist3(birth_g1)
% return

% histogram(g2_time, 'BinWidth', 0.1)
% title('Frequency distribution of G2 duration in FUCCI cells')
% xlabel('G2 duration')
% ylabel('Count')
% return

idx = 1;
% for time=0:h:120
%     t(idx)=time;
%     g1_count = 0;
%     g2_count = 0;
%     for row=1:size(A,1)
%         if(time>g1_start(row) && time+h<g1_end(row))
%             g1_count = g1_count + 1;
%         end
%         if(time>g2_start(row) && time+h<g2_end(row))
%             g2_count = g2_count + 1;
%         end
%     end
%     g1(idx)=g1_count;
%     g2(idx)=g2_count;
%     idx=idx+1;
%     hold off;
%     plot(t,100*g1./(g1+g2), 'b');
%     hold on;
%     plot(t,100*g2./(g1+g2), 'r');
%     title('Percent of FUCCI cells in G2 phase');
%     xlabel('Time (hrs)');
%     ylabel('Percentage');
%     drawnow;
% end

for time=0:h:120
    t(idx)=time;
    g1_count_A = 0;
    g2_count_A = 0;
    g1_count_B = 0;
    g2_count_B = 0;
    for row_censored=1:size(g1_start_censored,2)
        if(time>g1_start_censored(row_censored) && time+h<g1_end_censored(row_censored))
            if(censored_idx(row_censored)==1)
                g1_count_A = g1_count_A + 1;
            else
                g1_count_B = g1_count_B + 1;
            end
        end
        if(time>g2_start_censored(row_censored) && time+h<g2_end_censored(row_censored))
            if(censored_idx(row_censored)==1)
                g2_count_A = g2_count_A + 1;
            else
                g2_count_B = g2_count_B + 1;
            end
        end
    end
    g1_A(idx)=g1_count_A;
    g2_A(idx)=g2_count_A;
    g1_B(idx)=g1_count_B;
    g2_B(idx)=g2_count_B;
    idx=idx+1;
    hold off;
    total_pop=(g1_A + g2_A + g1_B + g2_B);
    total_g1=(g1_A + g1_B);
    total_g2=(g2_A + g2_B);
    total_popA=(g1_A + g2_A);
    total_popB=(g1_B + g2_B);
    plot(t,100*g1_A./total_pop, '-b');
    hold on;
    plot(t,100*g2_A./total_pop, '-r');
    plot(t,100*g1_B./total_pop, '--b');
    plot(t,100*g2_B./total_pop, '--r');
    plot(t,100*total_g1./total_pop, '-*b');
    plot(t,100*total_g2./total_pop, '-*r');
    plot(t,100*total_popA./total_pop, '-g');
    plot(t,100*total_popB./total_pop, '--g');
    title('Percent of FUCCI cells in G2 phase');
    xlabel('Time (hrs)');
    ylabel('Percentage');
    legend('% G1 Population A','% G2 Population A','% G1 Population B','% G2 Population B','% G1 Total','% G2 Total','% Population A','% Population B')
    drawnow;
end




end