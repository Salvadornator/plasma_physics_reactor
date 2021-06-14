function psitbxput(ver,shot,redo)

% PSITBXPUT		Computes and stores data in the RESULTS tree
%   PSITBXPUT(VER,SHOT[,REDO]) User should not invoke this function. It is
%   called automatically when accessing an empty or not up to date node
%   \RESULTS::PSITBX:*

if nargin < 3, redo = 1; end

if shot == -1 || shot >= 100000, tree = 'PCS';    
else,                            tree = 'RESULTS'; end
if isempty(mdsopen(tree,shot)), return, end
if strcmp(tree,'PCS')
 psi = psitbxtcv(shot);
 mdssetdef('\MGAMS.PSITBX')
else
 if ~redo & ver == mdsdata('\PSITBX:AS:VERSION_NUM'), return, end
 psi = psitbxtcv(shot,(0:200)/20);
 mdssetdef('\PSITBX')
end
if isempty(psi), return, end
fsd = psitbxp2p(psi,'FS');
rzs = psitbxg2g(fsd.grid,'C',fsd);
fsg = psitbxfsg(fsd);

put(tree,'RMAG'  ,fsd.rmag    ,{psi.t}  			  ,'m')
put(tree,'ZMAG'  ,fsd.zmag    ,{psi.t}  			  ,'m')
put(tree,'PSIMAG',fsd.psimag  ,{psi.t}  			  ,'Vs')
put(tree,'AS'    ,fsd.x       ,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put(tree,'RS'    ,rzs.x{1}    ,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put(tree,'ZS'    ,rzs.x{2}    ,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put(tree,'GRHO'  ,fsg.grho.x  ,{fsd.grid.x{1}		   ,psi.t},'m^-1')
put(tree,'SURF'  ,fsg.surf.x  ,{fsd.grid.x{1}		   ,psi.t},'m^2')
put(tree,'AREA'  ,fsg.area.x  ,{fsd.grid.x{1}		   ,psi.t},'m^2')
put(tree,'DAREA' ,fsg.darea.x ,{fsd.grid.x{1}		   ,psi.t},'m^2')
put(tree,'VOL'   ,fsg.vol.x   ,{fsd.grid.x{1}		   ,psi.t},'m^3')
put(tree,'DVOL'  ,fsg.dvol.x  ,{fsd.grid.x{1}		   ,psi.t},'m^3')

mdsclose

function put(tree,node,x,dim,unit)
n = length(dim);
expr = ['BUILD_SIGNAL(BUILD_WITH_UNITS($1,$2),'];
if n > 1, expr = [expr sprintf(',$%d',3:n+1)]; end
expr = [expr sprintf(',BUILD_WITH_UNITS($%d,"s"))',n+2)];
if strcmp(tree,'RESULTS')
 mdsput([node ':FOO'],expr,'x',mdscvt(x,'f'),unit,dim{:})
 mdsput([node ':VERSION_NUM'],psitbxver,'f')
else
 mdsput(node,expr,'x',mdscvt(x,'f'),unit,dim{:})
end
