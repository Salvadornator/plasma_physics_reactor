function f = thetastar(psi,g,fsd)

% PSITBXPSI/THETASTAR	straight field line angle

% THTSTAR = THETASTAR(PSI)
%
% MODIFICATION HISTORY:
%   22-03-01  HR  introduce eps to avoid division by zero.
%	18-01-10  JR  allow other grid
%				  allow fsd as input
%	10-02-10  JR  allow g with single rho. NB: fsd must be given in that case. (fsd on single rho doesn't work). NB2: fsd must contain enough flux surfaces (at least up to rho=1) (which ones exactly???).

if nargin < 2 | isempty(g) | ~isa(g,'psitbxgrid'),
	g = psitbxgrid('Flux','Grid','Default');
end
if nargin < 3 | isempty( fsd ),
	fsd = psitbxp2p( psi, 'FS', g );
end
if ~strcmpi( g.storage, 'Grid' ),
	error( 'g must be of type grid to allow integration along theta' )
end


% make sure psi contains rmag and zmag
%psi = psitbxmag(psi);
% can use 01 or any other format, since f is normalized by the full integral
% (any value specific to a format is cancelled out).
psi = psitbxp2p(psi,'01');

t   = psi.psitbxfun.t;
nt  = max( length(t), 1 );
nth = length(g.x{2});
k0  = repmat(reshape([1:nt],[1,1,nt]),size(g));

% cylindrical grid (with Time-Points storage)
gc  = psitbxg2g(g,'C',fsd);

%radial and vertical derivatives of psi
dpsidr = psitbxf2f(psi,gc,[1,0]);
dpsidz = psitbxf2f(psi,gc,[0,1]);

%keep only data field and check dimensionality wrt gc

dpsidr = dpsidr.x;
dpsidz = dpsidz.x;

s_dp = size( squeeze( dpsidr ) );
s_gc = size( squeeze( gc.x{1} ) );

if ~isequal( s_dp, s_gc )
	
	%dimension problem. Might happen if only 1 flux surface is given, for one or more times. In that case, gc and dpsidr are respectively line and column vectors.
	
	dpsidr = shiftdim( dpsidr, -1 ); %to get a singleton along dim 1
	dpsidz = shiftdim( dpsidz, -1 );
	
	%check
	if ~isequal( size( dpsidr ), s_gc )
		error( 'Dimension problems' )
	end
	
end

% integration of (R^2*J(theta,psi))^-1 over theta
% NB (JR): to follow the thesis formula, we should have (rmag-R)*dpsidr - (Z-zmag)*dpsidz
% but here both sign are inverted. Since F is then normalized by the full integral, it
% boils down to the same value.
F = cumsum( ...
		psitbxfun( ( ( gc.x{1} - psi.rmag( k0 ) + eps ).^2 + ( gc.x{2} - psi.zmag( k0 ) + eps ).^2 ) ./ ...
    				( gc.x{1} .* ( ( gc.x{1} - psi.rmag( k0 ) + eps ) .* dpsidr ...
	                 + ( gc.x{2} - psi.zmag( k0 ) + eps ) .* dpsidz ) ), ...
	              g,t), ...
  	2);

% normalize from 0 to 2pi
% then matches the first value of theta* with the first value of theta
% This last step is necessary to be consistent with the integration lower boundary.
f = F.*(2*pi./repmat(F.x(:,nth,:),[1,nth,1])) + g.x{2}(1);

