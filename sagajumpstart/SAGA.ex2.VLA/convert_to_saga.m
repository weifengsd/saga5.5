%Convert_To_SAGA.m
%
%
% Matlab script to convert individual frequency sio CSDM's
% to Peter's SAGA input format.  This calls write_covmat.m.
%
% set the number of elements
nels = 21;

% set the frequencies
%f = [49 64 79 94 112 130 148 166 201 235 283 338 388 ];
f = [388];

[junk,nf]=size(f);

a = ones(nels,nels*nf);

%
% load the data in
%
for i = 1:nf
 index = (i-1)*nels+1;
 % under UNIX
 %filein = ['covmat.' int2str(f(i))]
 %a(1:nels,index:index+nels-1)=conj(loadsio(filein,'complex'));
 % under Linux
 filein = ['fftData/csdm.',int2str(f(i)),'Hz.sio']
 d=sioread(filein,1,0,[1:42]);                 % read in data file
 c=complex(d(:,[1:2:41]),d(:,[2:2:42])); % put data in complex format

 a(1:nels,index:index+nels-1)=c;
end

% define the array element depths
%  for horizontal array
%depth = 213.0 .* ones(1,nels);
%
%  for vertical array
%depth = [ 94.125 99.755 105.38 111.00 116.62 122.25 127.88 ...
%          139.12 144.74 150.38 155.99 161.62 167.26 172.88 ...
%          178.49 184.12 189.76 195.38 200.99 206.62 212.25 ]
%
%  for reversed vertical array ( so tilt is implemented correctly )
depth = [ 212.25 206.62 200.99 195.38 189.76 184.12 178.49 
          172.88 167.26 161.62 155.99 150.38 144.74 139.12 
          127.88 122.25 116.62 111.00 105.38 99.755 94.125 ]
%
% write the data out;
%
fidout =  fopen('cov_dpss.in','a');
for ifreq=1:length(f)
   fprintf(fidout,' estimated covariance matrices using dpss\n');
   fprintf(fidout,' %f     0.000 dB\n',f(ifreq));
   fprintf(fidout,' %d\n',nels);
   fprintf(fidout,' %f\n',depth);
   cols = [ (ifreq-1)*nels+1:ifreq*nels ];
   Cx = a(:,cols); % cut one cov-matrix out of the bunch
   for row=1:nels
      for col=1:nels
         fprintf(fidout,'%10d%10d (%E,%E) \n',row,col,real(Cx(row,col)),imag(Cx(row,col)));
      end
   end
end
fprintf(fidout,'!  \n');
status=fclose(fidout);
