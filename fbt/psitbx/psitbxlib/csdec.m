function [M,taug,b] = csdec(t,x,xx,ec,d)

% CSDEC		Cubic B-spline decomposition
%   [M,TAUG,B] = CSDEC(T,X[,EC]) returns M such that M*Y are the
%   coefficients of the cubic B-spline combination that best approximates
%   Y(X) for the given knot sequence T. If T and X coincide, the exact
%   interpolation is obtained.
%   CSDEC(T,X,XX[,EC][,D]) returns M such that M*Y is the value of the D-th
%   derivative of this combination at XX.
%   EC specifies the edge conditions:
% 	'a' not-a-knot	f''' continuous at t(2) and t(end-1) (default)
% 	'n' natural	f''(t(1)) = f''(t(end)) = 0
% 	'p' periodic	f'(t(1)) = f'(t(end)) and f''(t(1)) = f''(t(end))
%	 's' symetric	f'(t(1)) = f''(t(end)) = 0
%	 't' symetric	f'(t(1)) = f'''(t(1)) = 0 (ill conditioned not recomended)
%	 'u' symetric	f'(t(1)) = f'''(t(1)) = f''(t(end)) = 0
%	 'v' symetric	f'(t(1)) = f'''(t(1)) = f'''(t(end)) = 0
%
%   For two variable functions Z(X,Y), CSDEC must be called for X and Y:
%   CSDEC(TX,X) * Z * CSDEC(TY,Y)'

% arguments
if nargin == 4
	if ischar(ec), d = 0; else d = ec; ec = 'a'; end
elseif nargin == 3
	if ischar(xx), ec = xx; else ec = 'a'; d = 0; end
elseif nargin < 3
	ec = 'a';
end
nt = length(t);
if nt < 2, error('Give more than one knot.'), end

% augment knot sequence
dt = mean(diff(t));
taug = [(-3:-1)'*dt+t(1);t(:);t(nt)+(1:3)'*dt];

% compute the augmented base function set
b = bspbase(taug, 4, x)';

% set up constraint matrix for specified edge conditions: T*c = 0
if nt == 2 % the solution is a straight line: f'' and f''' 0 at some mid-point
 T = [bspsum(taug, eye(4), (t(1)+t(2))/2, 2) * dt^2; ...
      bspsum(taug, eye(4), (t(1)+t(2))/2, 3) * dt^3];
elseif nt == 3 % the solution is a parabola: f''' 0 at some mid-points
 T = bspsum(taug, eye(5), [t(1)+t(2) t(2)+t(3)]/2, 3) * dt^3; 
elseif ec == 'n' % f''(t(1)) = f''(t(nt)) = 0
 T = bspsum(taug, eye(nt+2), t([1 nt]), 2) * dt^2;
elseif ec == 'a' % f''' equal at some mid-points
 T = [1 -1 0 0; 0 0 1 -1] * bspsum(taug, eye(nt+2), ...
     [t(1)+t(2) t(2)+t(3) t(nt-2)+t(nt-1) t(nt-1)+t(nt)]/2, 3) * dt^3;
elseif ec == 'p' % f'(t(1)) = f'(t(nt)) and f''(t(1)) = f''(t(nt))
 T = [1 -1 0 0; 0 0 1 -1] * ...
     [bspsum(taug, eye(nt+2), [t(1) t(nt)], 1) * dt;
      bspsum(taug, eye(nt+2), [t(1) t(nt)], 2) * dt^2];
elseif ec == 's' % f'(t(1)) = f''(t(end)) = 0
 T = [bspsum(taug, eye(nt+2), t(1), 1) * dt; ...
      bspsum(taug, eye(nt+2), t(nt), 2) * dt^2];
elseif ec == 't' % f'(t(1)) = f'''(t(1)) = 0
 T = [bspsum(taug, eye(nt+2), t(1), 1) * dt; ...
      bspsum(taug, eye(nt+2), (t(1)+t(2))/2, 3) * dt^3];
elseif ec == 'u' % f'(t(1)) = f'''(t(1)) = f''(t(nt)) = 0
 T = [bspsum(taug, eye(nt+2), t(1), 1) * dt; ...
      bspsum(taug, eye(nt+2), (t(1)+t(2))/2, 3) * dt^3;
						bspsum(taug, eye(nt+2), t(nt), 2) * dt^2];
elseif ec == 'v' % f'(t(1)) = f'''(t(1)) = f'''(t(nt)) = 0
 T = [bspsum(taug, eye(nt+2), t(1), 1) * dt; ...
      bspsum(taug, eye(nt+2), [t(1)+t(2) t(nt-1)+t(nt)]/2, 3) * dt^3];
else
 error('Invalid edge condition.')
end

% constraint regression
[q,~] = qr(T'); q = q(:,size(T,1)+1:end);
M = q/(b*q);

if nargin >= 3 && ~ischar(xx)
 M = bspsum(taug,M,xx,d);
end
