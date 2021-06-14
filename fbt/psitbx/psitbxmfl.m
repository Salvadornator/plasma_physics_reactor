function mfl = psitbxmfl(psi,RBt,gb,gl,pt)

% PSITBXPSI/PSITBXMFL   magnetic field lines
% PSITBXMFL(PSI,RBT[,GB,GL,PT)
% PSI a PSITBXFUN object containing the poloidal flux
% RBT [a time vector containing] RBt
% GB  a PSITBXGRID object defining the grid for field steps and the
%     bounding rectangle or that rectangle [R1 R2 Z1 Z2] or the mesh axis
%     of that grid {R Z}
% GL  a PSITBXGRID (grid, points or time points) for field line "start"
% PT  a bounding polygon [PTR PTZ] (the tiles)

% Modification history:
%   2012-02-17  HR  Keep correct sign when de-normalising poloidal flux.

% optional grid conversion
if nargin < 3, gb = []; end, gb = gf(gb,psi);
if nargin < 4, gl = []; end, gl = gf(gl,psi);
if nargin < 5, pt = []; end
if size(pt,1) == 2, pt = pt'; ptx = pt(:,1); pty = pt(:,2); end

% De-normalise poloidal flux
psi = psitbxp2p(psi,'*0');

% interpolate on new grids
if isa(gb,'psitbxgrid'), pb = psitbxf2f(psi.psitbxfun,gb);
else                     pb = psi.psitbxfun; gb = pb.grid; end
if isa(gl,'psitbxgrid'), pl = psitbxf2f(psi,gl);
else                     pl = psi.psitbxfun; gl = pl.grid; end
t = pl.t; nt = length(t); nl = size(pl); ll = prod(nl)/nt;

switch gl.storage(1)
 case 'G'
  [zl,rl] = meshgrid(gl.x{2},gl.x{1});
  rl = rl(:); zl = zl(:);
 case 'P'
  rl = repmat(gl.x{1}(:),1,nt);
  zl = repmat(gl.x{2}(:),1,nt);
 case 'T'
  rl = reshape(gl.x{1}(:),ll,nt); 
  zl = reshape(gl.x{2}(:),ll,nt);
end
plx = reshape(pl.x,ll,nt);
if length(plx) == 1, plx = [plx;plx]; end % for contourc with 1 contour

rb = gb.x{1}; zb = gb.x{2};

for kt = 1:nt
 pbx = pb.x(:,:,kt);
 
 cc = contourc(rb,zb,pbx',plx(:,kt));
 % limit if tile contour given
 if ~isempty(pt), cc = polygc(ptx,pty,cc); end
 % remove doublons, get size of contours
 l = 1; nc = 0; mc = 0;
 while l <= size(cc,2)
  k = find(all(diff(cc(:,l+1:l+cc(2,l)),1,2) == 0,1));
  cc(:,l+k) = []; cc(2,l) = cc(2,l) - length(k);
  nc = max(nc,cc(2,l));
  l = l + cc(2,l) + 1;
  mc = mc + 1;
 end
 % extract contours
 [rc,zc] = deal(NaN(mc,nc));
 pc = NaN(mc,1);
 l = 1;
 for k = 1:mc
  pc(k) = cc(1,l);
  rc(k,1:cc(2,l)) = cc(1,l+1:l+cc(2,l));
  zc(k,1:cc(2,l)) = cc(2,l+1:l+cc(2,l));
  l = l + cc(2,l) + 1;
 end
 % select closest contours
 lc = ones(1,ll);
 [R,Z,pr,pz] = deal(NaN(nc,ll));
 for k = 1:ll;
  l = find(pc == plx(k,kt));
  if ~isempty(l)
   l = l(imin(min((rc(l,:)-rl(k,kt)).^2+(zc(l,:)-zl(k,kt)).^2,[],2)));
   lc(k) = imin((rc(l,:)-rl(k,kt)).^2+(zc(l,:)-zl(k,kt)).^2);
   R(:,k) = rc(l,:);   
   Z(:,k) = zc(l,:);
  end
 end
  
 % find flux derivatives on contour points
 k = find(~isnan(R));
 px = psitbxfun(pbx,gb);
 tmp = psitbxgrid('c','p',{R(k),Z(k)});
 pr(k) = double(psitbxf2f(px,tmp,[1,0]));
 pz(k) = double(psitbxf2f(px,tmp,[0,1]));

 % toroidal direction step
 dr = diff(R,1,1);
 dz = diff(Z,1,1);
 pr = (pr(1:end-1,:) + pr(2:end,:))/2;
 pz = (pz(1:end-1,:) + pz(2:end,:))/2;
 dy = 2*pi*RBt(min(kt,length(RBt)))*(dz.*sign(pr) - dr.*sign(pz)) ./ (abs(pr) + abs(pz));
 
 % some contours are drawn opposite to the magnetic field -> negative dy
 for k = find(dy(1,:) < 0)
  l = sum(~isnan(R(:,k)));
  R(1:l,k) = R(l:-1:1,k);
  Z(1:l,k) = Z(l:-1:1,k);
  dy(1:l-1,k) = -dy(l-1:-1:1,k);
  lc(k) = l + 1 - lc(k);
 end
 
 % integrate along Y and shift so that Y=0 is input points
 Y = cumsum([zeros(1,ll);dy],1);
 Y = Y - repmat(Y(sub2ind([nc,ll],lc,1:ll)),nc,1);

end

mfl = psitbxgrid('O','T',{reshape(R,[nc,nl]),reshape(Y,[nc,nl]),reshape(Z,[nc,nl])},t);

function g = gf(g,psi)
switch class(g)
 case 'double' % grid given as a rectangle
  if ~isempty(g)
   rg = psi.psitbxfun.grid.x{1};
   zg = psi.psitbxfun.grid.x{2};
   tmp = g;
   rg = [tmp(1) rg(rg > tmp(1) & rg < tmp(2)) tmp(2)];
   zg = [tmp(3) zg(zg > tmp(3) & zg < tmp(4)) tmp(4)];
   g = psitbxgrid('Cylindrical','Grid',{rg,zg});
  end
 case 'cell' % grid given as {R,Z}
  g = psitbxgrid('Cylindrical','Grid',{sort(g{1}),sort(g{2})});
 case 'psitbxgrid' % grid given as a grid
 otherwise, error('Invalid grid specification')
end

