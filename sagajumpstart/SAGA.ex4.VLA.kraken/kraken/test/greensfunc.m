function [G] = greensfunc(r,zs,zr,U,kh,z,r_off); 
%--------------------------------------------------------------
% r = source-to-receiver range [in m]
% zs = source depth [in m]
% zr = vector of receiver depths [in m] (number of receivers x 1)
% U = matrix of eigenfunctions (number of depths x number of modes)
% kh = vector of complex eigenvalues (number of modes x 1)
% z = vector of sample depths of modal eigenfunctions [in m]
% r_off = vector of range offsets to account for tilt [in m]
%	  (number of receivers x 1)
%--------------------------------------------------------------
j = sqrt(-1);
tol1 = 0.1;
tol2 = 0.01;
[m,nummodes]=size(U);
if m ~= length(z)
  error('Dimension of mode matrix not consistent with mode sample depths')
end

% Find the array index corresponding to the source depth.
if any(abs(z-zs)<=tol1)
  izs = find(abs(z-zs)<=tol1);  % if length(izs)>1 use first element
  izs = izs(1);
else
  error('Source depth not sampled by eigenfunctions');
end

% Generate output eigenfunctions corresponding to sensor locations
numrcvrs = length(zr);
izr = [];
for i = 1:numrcvrs
  if any(abs(z-zr(i))<tol2)
    im = find(abs(z-zr(i))<tol2);
    izr = [izr;im(1)];
  end
end
Ur = U(izr,:);	% modes sampled by receiver depths

% Compute the modal amplitude.
s = (sqrt(2*pi*ones(nummodes,1))./sqrt(r*real(kh))).*U(izs,:)'.*exp(-j*kh*r);
%s = (sqrt(2*pi*ones(nummodes,1))./sqrt(r*real(kh))).*U(izs,:)' ...
%    .*exp(j*(kh*r + pi/4));

% We must accomodate for source-to-receiver range errors due to 
% tilt.  Modify U_s to account for differential phase across the 
% array by pre-multiplying by a tilt correction matrix.  Note:
% This correction factor is NOT accounted for in the denominator
% of the phase factor in the modal amplitude vector since the 
% correction has less effect in the denominator where its square
% root is taken.
for M = 1:nummodes
  for N = 1:numrcvrs
    phase(M,N) = exp(-j*kh(M,1)*r_off(N,1));
%    phase(M,N) = exp(j*(kh(M,1)*r_off(N,1) + pi/4));
  end;
  Phase = diag(phase(M,:));
  Ur_new(:,M) = Phase*Ur(:,M);
end;

% Compute the pressure field.
G = Ur_new*s;

