path(path,'C:/Users/wfeng/Documents/GitHub/saga5.5/matlab') % enter your saga matlab path here
plotsaga
vert
print('-dpsc','SAGAscatter.ps')
cmin = -15
cmax = 0
figure;
contsaga
vertC
print('-dpsc','SAGAcontour.ps')