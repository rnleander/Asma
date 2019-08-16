clear;
clf;
close all;
startTime = tic;

h=.01;
numAncestors = 2^14 % slow
%numAncestors = 2^10 % medium
%numAncestors = 2^6 % fast
maxGenerations = 7 % maximum number of descendent generations (not including root ancestor)
plotDDT = 0; % 0 disable, 1 deferred plot, 2 realtime plot
noGrowth = 1; % if 1, each cell has only a single descendent
cutdown=25; % Reduce the resolution of the G1 and G2 time/age plots by this factor
smooth = 0; % Use gaussian convolutional smoothing
model = 0; % Use 0 for driftdiffusion and 1 for exponential
startingOffset=10000; % Used to pad the initial t<0 region to avoid pointer arithmetic problems

set(gcf,'WindowStyle','docked')

% Run the simulation for a number of initial ancestors cells and collect a
% list of those cells
tic
if(plotDDT~=0)
    set(figure,'WindowStyle','docked')
    title('Discrete Simulation of Drift Diffusion Process')
    legend('X_1','X_2')
    hold on;
end
%fprintf("Simulation  0.00 percent complete\n")

descendents = cell(1,numAncestors);
pool = gcp;
numWorkers = pool.NumWorkers;
%numWorkers = 1;
clear myProgress;
myProgress("Simulation %5.2f percent complete\n",numAncestors, numWorkers);
parfor k=1:numAncestors
    ancestor = experiment(10000,0,maxGenerations,h,plotDDT,noGrowth,model,startingOffset); 
    descendents{k} = flattenDescendents(ancestor);
    myProgress("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bSimulation %5.2f percent complete\n",numAncestors, numWorkers);
end
if(plotDDT~=0)
    drawnow;
    lines = findobj(gca,'Type','line');
    legend([lines(1) lines(2)], 'X_1', 'X_2');
    xlabel('Time (hrs)');
    ylabel('Protein level');
    ylim([-0.5,1.5]);
end
toc

% Find the maximum time, maximum age and plot the IMT distribution
tic
fprintf('\nGathering cells and computing statistics\n');
allCells={};
%for k=1:numAncestors
%    allCells = cat(2,allCells,descendents{k});
%end
allCells = cat(2,descendents{:});

numAllCells = size(allCells,2)
maxT = 0;
maxAge = 0;
maxG1duration = 0;
maxG2duration = 0;
imt=[];
fixedIdx=0;
for cellIdx=1:numAllCells
    thisCell = allCells{cellIdx};
    if(thisCell.end > maxT)
        maxT = thisCell.end;
    end
    age = thisCell.end - thisCell.begin;
    if(age > maxAge)
        maxAge = age;
    end
    G1duration = thisCell.restrictionPoint - thisCell.begin;
    if(G1duration > maxG1duration)
        maxG1duration = G1duration;
    end
    G2duration = thisCell.end - thisCell.restrictionPoint;
    if(G2duration > maxG2duration)
        maxG2duration = G2duration;
    end
    imt(cellIdx) = thisCell.imt;
    if(thisCell.generation>0)
        fixedIdx = fixedIdx+1;
        imt_fixed(fixedIdx) = thisCell.imt;
    end
    
end
imt=imt*h;
set(figure,'WindowStyle','docked')
histogram(imt)
title('IMT Distribution including initial population');
xlabel('Time (hrs)');
ylabel('Density');
drawnow;
imt_fixed=imt_fixed*h;
set(figure,'WindowStyle','docked')
histogram(imt_fixed)
title('IMT Distribution not including initial population');
xlabel('Time (hrs)');
ylabel('Density');
drawnow;
maxT
maxAge
maxG1duration
maxG2duration
toc

% Plot the G1 and G2 percentages as a function of time and age
tic
tRange = floor(maxT/cutdown)
g1AgeRange = floor(maxG1duration/cutdown)
g2AgeRange = floor(maxG2duration/cutdown)
G1timeAge = zeros(tRange+1,g1AgeRange+1);
G1timeAgeNorm = zeros(tRange+1,g1AgeRange+1);
G2timeAge = zeros(tRange+1,g2AgeRange+1);
G2timeAgeNorm = zeros(tRange+1,g2AgeRange+1);
chkStep=round(numAllCells/50);
fprintf("Surface  0.00 percent complete\n");
% Accumulate krockner deltas for each cell life history for G1 and G2
% note - this uses the cutdown variable to reduce the resolution of the
% surface plot
for cellIdx=1:numAllCells
    thisCell = allCells{cellIdx};
    cutdownBeginning = floor(thisCell.begin/cutdown);
    for time=floor(thisCell.begin/cutdown)+1:floor(thisCell.restrictionPoint/cutdown)
        % HEY THIS NEEDS TO BE DOUBLECKECKED!!!
%          if(time<1)
%              time=1
%          end
       g1age = time - cutdownBeginning;
       G1timeAge(time, g1age) = G1timeAge(time, g1age) + 1;
    end
    for time=floor(thisCell.restrictionPoint/cutdown)+1:floor(thisCell.end/cutdown)
%          if(time<1)
%              time=1
%          end
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
if(smooth == 1)
    for i=1:1
        G1timeAgeNorm = imgaussfilt(G1timeAgeNorm,2);
        G2timeAgeNorm = imgaussfilt(G2timeAgeNorm,2);
    end
end

% Plot the initial conditions
set(figure,'WindowStyle','docked');
plot(G1timeAgeNorm(startingOffset/cutdown,:));
hold on;
plot(G1timeAgeNorm(startingOffset/cutdown,:));
title('G1 and G2 fractions at time t=0 (initial conditions)');
xlabel('time (hrs)');
ylabel('Density');
legend('G1','G2');

% Make plots for G1 and G2
[X,Y] = meshgrid(0:h*cutdown:maxT*h,0:h*cutdown:maxG1duration*h);
set(figure,'WindowStyle','docked');
surf(X,Y,G1timeAgeNorm','EdgeColor','none');
shading interp;
title('G1 fraction as a function of time and age');
xlabel('time (hrs)');
%xlims = xlim;
%xlim([0,xlims(2)]);
ylabel('age (hrs)');
zlabel('G1 Norm');
drawnow;
[X,Y] = meshgrid(0:h*cutdown:maxT*h,0:h*cutdown:maxG2duration*h);
set(figure,'WindowStyle','docked');
surf(X,Y,G2timeAgeNorm','EdgeColor','none');
shading interp;
title('G2 fraction as a function of time and age');
xlabel('time (hrs)');
%xlims = xlim;
%xlim([0,xlims(2)]);
ylabel('age (hrs)');
zlabel('G2 Norm');
drawnow;
toc

% Plot the G1 and G2 fractions as a function of time
tic
set(figure,'WindowStyle','docked')
g1 = zeros(1,maxT);
g2 = zeros(1,maxT);
chkStep=round(numAllCells/50);
fprintf("\nCounting  0.00 percent complete\n");
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
plot(0:h:(maxT-1)*h, g1./alive,'r');
hold on;
plot(0:h:(maxT-1)*h, g2./alive,'b');
%xlims = xlim;
%xlim([startingOffset*h,xlims(2)]);
drawnow;
toc

fprintf('\nTotal '); toc(startTime);

function myProgress(fmtstring,total,numWorkers)
    persistent n
    if isempty(n)
        n = 0;
    end
    fprintf(fmtstring, numWorkers*100*n/total);
    n = n+1;
end