function [ pltitl, freq, ck, z, phi, cpt, cst, rhot, deptht, cpb, csb, rhob, depthb, Nmedia, depth, rho ] =  ...
  read_modes( filename, modes );

% read_modes
% useage:
%    [ pltitl, freq, ck, z, phi ] = read_modes( filename, modes )
% read the modes produced by KRAKEN
% filename should include the extension
% modes is an optional vector of mode indices

% mbp, May 2001

%determine type of file:

switch filename
case 'MODFIL'
    extension = '.mod';
otherwise
    endchar = length( filename );
    extension = lower( filename( endchar - 3 : endchar ) );
end

switch extension
case '.mod' % is in binary format
  if nargin == 1
     [ pltitl, freq, ck, z, phi, cpt, cst, rhot, deptht, cpb, csb, rhob, depthb, Nmedia, depth, rho ] = read_modes_bin( filename );
  else
     [ pltitl, freq, ck, z, phi, cpt, cst, rhot, deptht, cpb, csb, rhob, depthb, Nmedia, depth, rho ] = read_modes_bin( filename, modes );
  end
case '.moa'
  if nargin == 1
    [ pltitl, freq, ck, z, phi, cpt, cst, rhot, deptht, cpb, csb, rhob, depthb, Nmedia, depth, rho ] = read_modes_asc( filename );
  else
    [ pltitl, freq, ck, z, phi, cpt, cst, rhot, deptht, cpb, csb, rhob, depthb, Nmedia, depth, rho ] = read_modes_asc( filename, modes );
  end
otherwise
   error( 'unrecognized file extension' )
end
