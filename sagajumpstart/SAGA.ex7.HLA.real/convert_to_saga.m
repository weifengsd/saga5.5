%Convert_To_SAGA.m
%
%
% Matlab script to convert individual frequency sio CSDM's
% to Peter's SAGA input format.  This calls write_covmat.m.
%
% set the number of elements
nels = 27;

% set the frequencies
f = [49 64 79 94 112 130 148 166 201 235 283 338 388];

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
% filein = ['fftData/csdm.',int2str(f(i)),'Hz.sio']
 filein = ['csdm/csdm.',int2str(f(i)),'Hz.sio']
 d=sioread(filein,1,0,[1:54]);                 % read in data file
 % input data needs to be un-conj ( was complex conj for NBF conformability )
 % c=conj(complex(d(:,[1:2:53]),d(:,[2:2:54]))); % put data in complex format, conj
 c=complex(d(:,[1:2:53]),d(:,[2:2:54])); % put data in complex format
 a(1:nels,index:index+nels-1)=c;
end

% define the array element depths
depth = 213.0 .* ones(1,nels);
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
