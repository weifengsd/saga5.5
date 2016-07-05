function status=write_covmat(freqvec,Cov)
%write_cov_dpss(freqvec,Cov)
%
% Write Covariance Matrices to file cov_dpss.in
%
% freqvec = vector of frequencies in Hertz
% Cov     = matrix of concatenated covariance matrices
%           this must be of size [27, length(freqvec)*27]
[n,m]= size(Cov);
if(m ~= n*length(freqvec))
   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
           size(Cov),length(freqvec));
   error('sizes of Cov and freqvec are not compatible');
end
%depth=linspace(18.72,112.72,n);
n
depth = 213.0 .* ones(1,n);

fidout =  fopen('cov_dpss.in','a');

for ifreq=1:length(freqvec)
   fprintf(fidout,' estimated covariance matrices using dpss\n');
   fprintf(fidout,' %f     0.000 dB\n',freqvec(ifreq));
   fprintf(fidout,' %d\n',n);
   fprintf(fidout,' %f\n',depth);
   cols = [ (ifreq-1)*n+1:ifreq*n ];
   Cx = Cov(:,cols); % cut one cov-matrix out of the bunch
   for row=1:n
      for col=1:n
         fprintf(fidout,'%10d%10d (%E,%E) \n',row,col, ...
                 real(Cx(row,col)),imag(Cx(row,col)));
      end
   end
end
fprintf(fidout,'!  \n');
status=fclose(fidout);
unix('ls -l cov_dpss.in');

