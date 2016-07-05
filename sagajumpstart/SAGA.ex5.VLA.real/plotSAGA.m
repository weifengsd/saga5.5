clear; close all;

% path(path,'/kingcone0/dave/SAGA/SwellEx96/saga51/matlab') % enter your saga matlab path here
plotsaga; % Input: vert
print('-dpsc','SAGAscatter.ps')
cmin = -15
cmax = 0
figure;
contsaga; % Input: vertC
print('-dpsc','SAGAcontour.ps')
figure;
contsaga; % Input: vertC2
print('-dpsc','SAGAcontour2.ps')
