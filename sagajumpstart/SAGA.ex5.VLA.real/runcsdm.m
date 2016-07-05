%%  script to generate covariance matrix for SAGA with real data using a VLA 
clear;
prog_dir = [pwd '/fftData'];
path(prog_dir,path);
% remove old data, if any
if exist('cov.in','file')
    delete (which('cov.in'));
end

%% now compute the covariance matrix from data and convert covariance matrix to SAGA format 
% cd fftData
compute_csdm; % program file name to compute and write out covariance matrix
% cd ..


