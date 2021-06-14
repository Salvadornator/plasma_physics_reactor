function psi = psitbxmag(psi)

% PSITBXPSI/PSITBXMAG

% Modif log:
%
%   12.07.10  JR  Make sure that a (dumb) call with '01' input keeps the
%                 sign of psimag if .psimag filled. This case should not
%                 happen, but psitbxmag is now more robust (in case
%                 psitbxpsi is modified later on).

if ~isempty(psi.rmag) & strcmp(psi.format,'01'), return, end

x = psi.psitbxfun.x;
g = psi.psitbxfun.grid;
blnNegpm = false; %JR 12.07.10
switch(psi.format)
 case '-0', 
     x = -x;
     blnNegpm = true; %JR 12.07.10
 case '01', 
     x = 1 - x;
     if ~isempty( psi.psimag ) & mean( psi.psimag ) < 0
         blnNegpm = true; %JR 12.07.10
     end
 case 'FS', error('"Flux-Surfaces" can not be normalised')
end

% alloc
[nr,nz,nt] = size(x);
[psi.rmag,psi.zmag,psi.psimag] = deal(repmat(NaN,1,nt));

Mx = g.csp.Mx;
xk = g.csp.xk;
for kt = 1:nt

 % spline coefficients
 xx = x(:,:,kt);
 c = Mx{1}*xx*Mx{2}';
 
 % initial guess, but not on the grid edges
 kr0 = 1; kr1 = nr; kz0 = 1; kz1 = nz;
 while 1
  [kr,kz] = find(xx == max(xx(:))); kr = kr(1); kz = kz(1); % in case of two maxima
  if     kz == kz0, kz0 = kz0 + 1; xx(:,kz) = NaN;
  elseif kz == kz1, kz1 = kz1 - 1; xx(:,kz) = NaN;
  elseif kr == kr0, kr0 = kr0 + 1; xx(kr,:) = NaN;
  elseif kr == kr1, kr1 = kr1 - 1; xx(kr,:) = NaN;
  else, break
  end
 end
 r = g.x{1}(kr); z = g.x{2}(kz);

 % Gauss-Newton iteration
 for kiter = 1:psi.iter
 
  cr0 =  bspsum(xk{1},c,  r,0,0,1)';
  cr1 =  bspsum(xk{1},c,  r,1,0,1)';
  dp  = [bspsum(xk{2},cr1,z,0,0,1);
         bspsum(xk{2},cr0,z,1,0,1)];
  d2p = [bspsum(xk{2},bspsum(xk{1},c,r,2)',z,0,0,1), ...
         bspsum(xk{2},cr1,z,1,1); ...
         0, ...
	 bspsum(xk{2},cr0,z,2,0,1)];
  d2p(2,1) = d2p(1,2);
  drz = d2p \ dp;
  r = r - drz(1); z = z - drz(2);
  if max(abs(drz)) < psi.tol, break, end
 end
 
 if kiter < psi.iter
  psi.rmag(kt) = r; psi.zmag(kt) = z;
  psi.psimag(kt) = bspsum(xk{2},bspsum(xk{1},c,r,0,0,1)',z,0,1,1);
  x(:,:,kt) = 1 - x(:,:,kt) / psi.psimag(kt);
 else
  x(:,:,kt) = NaN;
 end
end

if blnNegpm, psi.psimag = -psi.psimag; end %JR 12.07.10. Now uses blnNegpm instead of test on input format.
psi.format = '01';
psi.psitbxfun = psitbxfun(x,g,psi.psitbxfun.t);
