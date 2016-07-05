clear;
path(path,'$SAGADIR/matlab') % enter your saga matlab path here

plotsaga % vert
print('-dpsc','SAGAscatter.ps')
cmin = -10
cmax =  -0
figure;
contsaga % vertC

print('-dpsc','SAGAcontour.ps')
figure;
contsaga % vertC2

print('-dpsc','SAGAcontour2.ps')

% # clean up
% rm fort.* snap*dat cov.in
% rm vert*.mat vert*bdr vert*cdr vert*.m vert*.out vert.cor vert.ext
% rm vert.p* vert.trf results.m *.obs