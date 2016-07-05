
% function to readin data and calculate CSDM matrix K for each frequency of
% interest

clear all

T0 = 0.1; % Start at T0 minute.
% Ntotal = 22*2^13;
Ntotal = 10*2^13;
N_elt = 21;
fs = 1500; 
Nfft = 2^14; 
segment_length = 2^13;
f_bin = fs/Nfft;
frange = 0:f_bin:0.5*Nfft*f_bin;
delta_t = segment_length/fs/60;
start_t = 0;
end_t = start_t + Ntotal/fs/60;
trange = start_t:delta_t:end_t;
data = sioread('J1312355.vla.21els.sio', T0*60*fs, Ntotal, 1:N_elt); % Swell-Ex96 data
depth = [212.25 206.62 200.99 195.38 189.76 184.12 178.49 172.88 167.26 161.62 155.99 150.38 ...
         144.74 139.12 127.88 122.25 116.62 111.00 105.38 99.755 94.125]; % depth of each hydrophone 

% Perform fft with 50% overlap using Kaiser-Bessel window (beta = 7.85) for
% each channel

beta = 7.85;
window = kaiser(segment_length,beta);
Noverlap = 0.5*segment_length; % specify 50% overlap 
N_chunks = 2*ceil(length(data)/segment_length)-1; % no. of segments with 50% overlap

% window each chunk and perform FFT with 50% overlap
for i_elt = 1:1:N_elt
    for i_chunk = 0:1:N_chunks-1
        start_index = i_chunk*0.5*segment_length+1;
        end_index = start_index+segment_length-1;
        chunk = data(start_index:end_index,i_elt).*window;
        FFTdata(i_elt,:,i_chunk+1) = fft(chunk,Nfft); 
    end
end    

FFTdata = FFTdata(:,1:Nfft/2,:);
FFTdata_mean = mean(FFTdata,1);
% figure;
% imagesc(trange,frange,20*log10(abs(squeeze(FFTdata_mean))));
% xlabel('time (min)');
% ylabel('Freq (Hz)');
% axis([start_t end_t 0 450]);
% caxis([-85 0]);

% calculate frequency bin of interest
f_interest = [49 64 79 94 112 130 148 166 201 235 283 338 388]; 
f_index = round(f_interest./f_bin); % find index of f_interest along frequency axis

% Search FT_data to extract 21 frequency bins at each frequency of interest (+- 10 bins around center frequency) 
% to take into account Doppler. 
% Calculate the power vs. time for each bin and separate each frequency set of bins into a file containing the 
% peak bin for each time snapshot.

FT_dataset = zeros(N_elt,length(f_index),N_chunks);
dk = zeros(N_elt,N_chunks,length(f_index));

for i_elt = 1:1:N_elt
    for i_freq = 1:1:length(f_index)
        for i_chunk = 1:1:N_chunks
            FT_dataset(i_elt,i_freq,i_chunk) = max(FFTdata(i_elt,f_index(i_freq)-10:1:f_index(i_freq)+10,i_chunk));
            dk(i_elt,i_chunk,i_freq) = FT_dataset(i_elt,i_freq,i_chunk);
        end
    end    
end

% now FT_dataset contains the set with max value for each elt, at the
% 13 frequency bins of interest, and time chunk (21*13*43).

K = zeros(N_elt,N_elt,length(f_index));    

for i_freq = 1:1:length(f_index)
    for i_chunk = 1:1:N_chunks
        K(:,:,i_freq) = K(:,:,i_freq) + dk(:,i_chunk,i_freq)*(dk(:,i_chunk,i_freq))';
    end
    K(:,:,i_freq) = K(:,:,i_freq)/N_chunks;
end

cov_file = write_covmat(f_interest,K,depth); % write file to cov.in in saga format.
eval(['!/bin/mv cov.in ../cov.in']); % move cov.in into the directory above where the .dat file is
