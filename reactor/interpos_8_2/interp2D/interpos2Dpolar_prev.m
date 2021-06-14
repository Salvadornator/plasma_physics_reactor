function [f_out,dfdrho_out,dfdtheta_out,d2fdrho2_out,d2fdtheta2_out,d2fdrhodtheta_out,varargout]= ...
    interpos2Dpolar(rho_in,theta_in,f_in,rho_out,theta_out,tension_rho,tension_theta,nbc_rho,kextrapol,varargin);
%
% [f_out,dfdrho_out,dfdtheta_out,d2fdrho2_out,d2fdtheta2_out,d2fdrhodtheta_out,varargout]= ...
%        interpos2Dpolar(rho_in,theta_in,f_in,rho_out,theta_out,tension_rho,tension_theta,kextrapol,varargin);
%
%
% Assume 1-D input meshes rho(i), theta(j) and 2-D input function f_in(i,j)
%
% Compute the values of f at all points (rho_out, theta_out) which can be a scalar, 1D or 2D arrays.
%
% It can compute also the first or second derivatives
%
% The interpolation is performed with two successive series of 1-D cubic spline interpolation with tension.
% First one calculates the values at (rho(i),theta_out), using periodic B.C and then one interpolates at the rho_out values.
%
% For the rho interpolation, boundary conditions and extrapolation rule are important. If the value at rho=0 is given, it can be used (with a non-zero tension). 
% This is the best, except if zero radial derivative is expected. The B.C. follow the interpos rules:
%        nbc_rho=[2 2] (and the routine uses f(1,:) and f(end,:)): Values at both ends to be used
%        nbc_rho=[1 2] df/drho(rho=0)=0, f(end,:) used
%        nbc_rho=[0 0] standard B.C for splines with zero second derivatives at each end.
%
%     extrapolate according to kextrapol (only linear implemented at this stage):  
%     kextrapol = 0 => constant
%     kextrapol = 1 => linear
%     kextrapol = 2 => quadratic
%     kextrapol = 3 => cubic (only safe if does not extrapolate very far)
%
f_out = [];
dfdrho_out = [];
dfdtheta_out = [];
d2fdrho2_out = [];
d2fdtheta2_out = [];
d2fdrhodtheta_out = [];
%
% Check inputs:
%
% 
if prod(size(rho_in)) ~= length(rho_in)
  warning('expects 1-D rho_in array')
  return
end
if prod(size(theta_in)) ~= length(theta_in)
  warning('expects 1-D rho_in array')
  return
end
if numel(f_in) ~= numel(rho_in).*numel(theta_in)
  disp('bad number of input elements, check input arrays')
  return
end
if numel(rho_out) ~= numel(theta_out)
  disp('bad number of out mesh points, check rho_out, theta_out arrays')
  return
end
if length(nbc_rho) ~= 2
  disp('bad number boundary conditions for rho, needs 2, check nbc_rho')
  return
end
% Check if theta end points doubled
theta_in_eff = theta_in;
f_in_eff = f_in;
if abs(theta_in_eff(end)-theta_in_eff(1)-2.*pi) <=1e-4*min(diff(theta_in_eff))
  % disp('Seems 2pi theta point  given, do not use it')
  theta_in_eff = theta_in_eff(1:end-1);
  f_in_eff = f_in_eff(:,1:end-1);
end
rho_in_eff = rho_in;
nb_in = length(rho_in_eff);

% Needs 2-D rho_out, theta_out for generality
if numel(rho_out) == max(size(rho_out))
  ndim1=length(rho_out);
  ndim2=1;
  rho_out_eff = reshape(rho_out,ndim1,ndim2);
  theta_out_eff = reshape(theta_out,ndim1,ndim2);
else
  ndim1 = size(rho_out,1);
  ndim2 = size(rho_out,2);
  rho_out_eff = rho_out;
  theta_out_eff = theta_out;
end
f_out = NaN * ones(ndim1,ndim2);
dfdrho_out = NaN * ones(ndim1,ndim2);
dfdtheta_out = NaN * ones(ndim1,ndim2);
d2fdrho2_out = NaN * ones(ndim1,ndim2);
d2fdtheta2_out = NaN * ones(ndim1,ndim2);
d2fdrhodtheta_out = NaN * ones(ndim1,ndim2);
keyboard
for j=1:ndim2
  for i_rho_in=1:nb_in
    f_rhoin_thetaout(i_rho_in,:) = interpos(theta_in_eff-pi,f_in_eff(i_rho_in,:),theta_out_eff(:,j),tension_theta,[-1 ],[2.*pi]);
  end
  for i=1:ndim1
    if nbc_rho(1) <= 1
      ybc_rho(1) = 0
    elseif nbc_rho(1) == 2
      ybc_rho(1) = f_rhoin_thetaout(1,i);
    end
    if nbc_rho(2) <= 1
      ybc_rho(1) = 0
    elseif nbc_rho(2) == 2
      ybc_rho(2) = f_rhoin_thetaout(end,i);
    end
    f_out_eff(i,j) = interpos(rho_in_eff,f_rhoin_thetaout(:,i),rho_out_eff(i,j),tension_rho,nbc_rho,ybc_rho);
  end
end
  
