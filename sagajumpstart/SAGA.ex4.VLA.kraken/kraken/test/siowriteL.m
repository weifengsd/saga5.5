function status=siowriteL(filename,x);
%-----------------------------------------------------------------------
% siowrite
%
% This program runs under windows, unix, and macs -- but
% writes "Unix" sio files.
%
% function status=siowrite(filename,array);
%
% filename = sio file to create
% array = variable to write to sio file
%
% written by James Murray (Final Version: July 2001)
% modified by David Ensberg to speed up code (May 2003)
%
%-----------------------------------------------------------------------

% For Now, fix the record length at 8192
rl = 8192;
% also, output floats only
sl = 4;
type = 'float';

tp = size(x);			
np = tp(1);			% number of points (per channel)
nc = tp(2);			% number of channels
nr = ceil(np/rl*sl)*nc;		% total number of records (not including header)
nrf = full(np/rl*sl)*nc;        % total number of full records(excluding header)
ppr = rl/sl;			% number of points per record

%endian = 'b';
endian = 'l';
fid=fopen(filename,'w',endian);

% Write Header Record
a = zeros(2048,1);
a(1) = 10;		% id value
a(2) = nr;		% number of records
a(3) = rl;		% record length in bytes
a(4) = nc;		% number of channels
a(5) = sl;		% bytes/point
a(6) = 1;		% real output
a(7) = np;		% number of points
a(8) = 32677;		% check point

% write the filename as their ASCII values... I know I could do this better
for jj = 1:23
  s(jj) = ' ';
end
for jj = 1:min(24,length(filename))
 s(jj) = filename(jj);
end
 a(9)=hex2dec([dec2hex(s(1)) dec2hex(s(2)) dec2hex(s(3)) dec2hex(s(4))]);
 a(10)=hex2dec([dec2hex(s(5)) dec2hex(s(6)) dec2hex(s(7)) dec2hex(s(8))]);
 a(11)=hex2dec([dec2hex(s(9)) dec2hex(s(10)) dec2hex(s(11)) dec2hex(s(12))]);
 a(12)=hex2dec([dec2hex(s(13)) dec2hex(s(14)) dec2hex(s(15)) dec2hex(s(16))]);
 a(13)=hex2dec([dec2hex(s(17)) dec2hex(s(18)) dec2hex(s(19)) dec2hex(s(20))]);
 a(14)=hex2dec([dec2hex(s(21)) dec2hex(s(22)) dec2hex(s(23)) dec2hex(s(23))]);


status = fwrite(fid,a,'int32');
  if status == -1
    error(ferror(fid))
  end

% write out full records 
nr1 = nrf/nc;
count = 1;
tmp=zeros(ppr,1);
for jj = 1:nr1
  countEnd = count+ppr-1;
  for ii = 1:nc
    status = fwrite(fid,x(count:countEnd,ii),'real*4');
  end
  count = count + ppr;
end

% now write out partial records ( if any )
if   nr ~= nrf ,
  countEnd = count+min(ppr,np-count);
  for ii = 1:nc
    t= x(count:countEnd,ii);
    tmp = [t',zeros(ppr-length(t),1)']'; % pad out to length of full record
    status = fwrite(fid,tmp,'real*4');
  end
  count = count + ppr;
end
   

fclose(fid);
