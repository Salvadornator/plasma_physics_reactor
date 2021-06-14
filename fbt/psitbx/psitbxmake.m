function psitbxmake(mode)

% PSITBXMAKE   Make file for PSITBX
%   psitbxmake psitbxmex   build PSITBX specific mex's
%   psitbxmake psitbxlib   build other mex's in psitbxlib directory
%   psitbxmake distrib     create links and /tmp/psitbx.tar for distribution
%   psitbxmake install     to be run for stand alone installation

if nargin == 0, mode = 'psitbxmex'; end

libfiles = { 'bspbase.m' 'bspbase.c' 'bspbasemex.c' 'bspsum.m' 'bspsummex.c' 'csdec.m' 'xmtimes.m' 'xop.m' };

switch computer
 case 'IBM_RS'
  mexoptions = { 'COPTIMFLAGS=-qsmp -qarch=auto -qtune=auto -DNDEBUG' 'LDOPTIMFLAGS=-s -qsmp' };
 otherwise
  mexoptions = {}; 
end

switch lower(mode)
 
 case 'psitbxmex'
  mex(mexoptions{:},'-outdir','@psitbxdcd/private/','@psitbxdcd/private/intersect2dmex.c')
  
 case 'psitbxlib'
  mex(mexoptions{:},'-outdir','psitbxlib','psitbxlib/bspbasemex.c','psitbxlib/bspbase.c')
  mex(mexoptions{:},'-outdir','psitbxlib','-DBSPPAR','-DBSPSUM4','psitbxlib/bspsummex.c')
  
 case 'distrib'
  for s = libfiles
   t = which(s{1});
   if ~isempty(t), unix(['ln -sf ' t ' psitbxlib']); end
  end
  unix('ln -sf ../html');
  unix('tar -hcf /tmp/psitbx.tar .');
  
 case 'install'
  psitbxmake psitbxmex
  psitbxmake psitbxlib
  copyfile psitbxlib/* .
  
 otherwise
  error('Invalid MODE')
end
