clear;
clf;

h=.01;
%numAncestors = 2^12 % slow
numAncestors = 2^8 % fast
maxGenerations = 7 % maximum number of descendent generations (not including root ancestor)
plotDDT = 0; % 0 disable, 1 deferred plot, 2 realtime plot
noGrowth = 1; % if 1, each cell has only a single descendent
cutdown=25; % Reduce the resolution of the G1 and G2 time/age plots by this factor


% Run the simulation for a number of initial ancestors cells and collect a
% list of those cells
if(plotDDT~=0)
    set(figure,'WindowStyle','docked')
    hold on;
end
fprintf("Simulation  0.00 percent complete\n")
descendents = cell(1,numAncestors);
for k=1:numAncestors
    ancestor = experiment(0,0,maxGenerations,h,plotDDT,noGrowth); 
    descendents{k} = flattenDescendents(ancestor);
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    fprintf("%5.2f percent complete\n", 100*k/numAncestors);
end
if(plotDDT==1)
    drawnow;
end
allCells={};
for k=1:numAncestors
    allCells = cat(2,allCells,descendents{k});
end
numAllCells = size(allCells,2)


% Find the maximum time, maximum age and plot the IMT distribution
maxT = 0;
maxAge = 0;
imt=[];
for cellIdx=1:numAllCells
    thisCell = allCells{cellIdx};
    if(thisCell.end > maxT)
        maxT = thisCell.end;
    end
    age = thisCell.end - thisCell.begin;
    if(age > maxAge)
        maxAge = age;
    end
    imt(cellIdx) = thisCell.imt;
end
imt=imt*h;
set(figure,'WindowStyle','docked')
histogram(imt)
title('IMT Distribution');
xlabel('Time (hrs)');
ylabel('Density');
drawnow;


% Plot the G1 and G2 percentages as a function of time and age
G1timeAge = zeros(floor(maxT/cutdown)+1,floor(maxAge/cutdown)+1);
G1timeAgeNorm = zeros(floor(maxT/cutdown)+1,floor(maxAge/cutdown)+1);
G2timeAge = zeros(floor(maxT/cutdown)+1,floor(maxAge/cutdown)+1);
G2timeAgeNorm = zeros(floor(maxT/cutdown)+1,floor(maxAge/cutdown)+1);
chkStep=round(numAllCells/50);
fprintf("Surface  0.00 percent complete\n");
% Accumulate krockner deltas for each cell life history for G1 and G2
% note - this uses the cutdown variable to reduce the resolution of the
% surface plot
for cellIdx=1:numAllCells
    thisCell = allCells{cellIdx};
    cutdownBeginning = floor(thisCell.begin/cutdown);
    for time=floor(thisCell.begin/cutdown)+1:floor(thisCell.restrictionPoint/cutdown)
       g1age = time - cutdownBeginning;
       G1timeAge(time, g1age) = G1timeAge(time, g1age) + 1;
    end
    for time=floor(thisCell.restrictionPoint/cutdown)+1:floor(thisCell.end/cutdown)
       g2age = time - floor(thisCell.restrictionPoint/cutdown);
       G2timeAge(time, g2age) = G2timeAge(time, g2age) + 1;
    end
    if(0==mod(cellIdx,chkStep))
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        fprintf("%5.2f percent complete\n", floor(100*cellIdx/numAllCells));
    end
end
% Normalize
for t=1:floor(maxT/cutdown)
    %totalTimeAgeIntegral = sum(G1timeAge(t,:)+G2timeAge(t,:));
    G1timeAgeNorm(t,:) = G1timeAge(t,:)./sum(G1timeAge(t,:));
    G2timeAgeNorm(t,:) = G2timeAge(t,:)./sum(G2timeAge(t,:));
end
% Gaussian blur to hide discrete noise
for i=1:1
    G1timeAgeNorm = imgaussfilt(G1timeAgeNorm,2);
    G2timeAgeNorm = imgaussfilt(G2timeAgeNorm,2);
end
% Make plots for G1 and G2
[X,Y] = meshgrid(0:h*cutdown:maxT*h,0:h*cutdown:maxAge*h);
set(figure,'WindowStyle','docked');
surf(X,Y,G1timeAgeNorm','EdgeColor','none');
shading interp;
title('G1 fraction as a function of time and age');
xlabel('time (hrs)');
ylabel('age (hrs)');
zlabel('G1 Norm');
drawnow;
set(figure,'WindowStyle','docked');
surf(X,Y,G2timeAgeNorm','EdgeColor','none');
shading interp;
title('G2 fraction as a function of time and age');
xlabel('time (hrs)');
ylabel('age (hrs)');
zlabel('G2 Norm');
drawnow;


% Plot the G1 and G2 fractions as a function of time
set(figure,'WindowStyle','docked')
g1 = zeros(1,maxT);
g2 = zeros(1,maxT);
chkStep=round(numAllCells/50);
fprintf("Counting  0.00 percent complete\n");
for cellIdx=1:numAllCells
    thisCell = allCells{cellIdx};
    for i=thisCell.begin+1:thisCell.restrictionPoint
        g1(i) = g1(i) + 1;
    end
    for i=thisCell.restrictionPoint+1:thisCell.end
        g2(i) = g2(i) + 1;
    end
    if(0==mod(cellIdx,chkStep))
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        fprintf("%5.2f percent complete\n", floor(100*cellIdx/numAllCells));
    end
%     alive = g1+g2;
%     hold off;
%     plot(0:h:(maxT-1)*h, g1./alive,'r')
%     hold on;
%     plot(0:h:(maxT-1)*h, g2./alive,'b')
%     drawnow;
end
alive = g1+g2;
plot(0:h:(maxT-1)*h, g1./alive,'r')
hold on;
plot(0:h:(maxT-1)*h, g2./alive,'b')
drawnow;



