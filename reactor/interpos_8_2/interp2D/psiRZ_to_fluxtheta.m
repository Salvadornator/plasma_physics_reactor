load equil.mat

%
% Use profiles_2d psi(grid) assuming grid is dim1=R, dim2=Z rectangular mesh
% Find flux surfaces inside the LCFS (at this stage but could interpolate outside)
% Flux surfaces are chosen from a rhopsi_norm_in mesh or a rhotor_norm_in (here using nflux points)
% The theta_mesh is equidistant using ntheta points
%
% Interpolates also br, bz and bphi on the flux surface mesh
%
% Needs also:
% plasma boundary (Rbnd, Zbnd) taken from equil.eqgeometry.boundary.r/z
% plasma axis from : equil.global_param.mag_axis.position.r/z
% Either: psi_axis/psi_edge from equil.global_param.psi_ax/bound to normalize equil.profiles_2d.psi
%         phi_edge from equil.profiles_1d.phi(end) to normalize equil.profiles_2d.phi
%
RR=equil.profiles_2d.grid.dim1;
ZZ=equil.profiles_2d.grid.dim2;
Rbnd=equil.eqgeometry.boundary.r;
Zbnd=equil.eqgeometry.boundary.z;

Raxis = equil.global_param.mag_axis.position.r;
Zaxis = equil.global_param.mag_axis.position.z;

nflux_out=80;
ntheta_out=121;

use_psi=1; % 1 to use psi, 0 to use phi
if use_psi
  psi_axis = equil.global_param.psi_ax;
  psi_edge = equil.global_param.psi_bound;
  psi_rz=equil.profiles_2d.psi;
  flux_rz_norm = (psi_rz-psi_axis) ./ (psi_edge-psi_axis);
  flux_norm_out=linspace(0.,1.,nflux_out); % insert here desired psi_norm_out mesh for flux surfaces
else
  phi_edge = equil.profiles_1d.phi(end);
  phi_rz = equil.profiles_2d.phi .* equil.global_param.toroid_field.r0.^2 .* equil.global_param.toroid_field.b0; % norm missing in CHEASE
  flux_rz_norm = phi_rz ./ phi_edge; % phi_axis=0 by definition
  flux_norm_out=linspace(0.,1.,nflux_out); % insert here desired phi_norm_out mesh for flux surfaces
end

% 2D quantity to get flux surface in (:,:,1) and then quantities to interpolate in farray_in(:,:,2:end)
farray_in(:,:,1) = flux_rz_norm;
farray_in(:,:,2) = equil.profiles_2d.br;
farray_in(:,:,3) = equil.profiles_2d.bz;
farray_in(:,:,4) = equil.profiles_2d.bphi;

% numerics choices
tension_default = -0.1;

% Compute a rho,theta mesh of plasma boundary. Make theta in 0,2pi and increasing for fitting
rhobnd=sqrt((Rbnd-Raxis).^2 + (Zbnd-Zaxis).^2);
thetabnd=atan2(Zbnd-Zaxis,Rbnd-Raxis);
ii=find(thetabnd<0);
thetabnd(ii) = thetabnd(ii) + 2*pi;
[thetabndsort,isort] = sort(thetabnd);
rhobndsort = rhobnd(isort);

% theta mesh
thetamesh=linspace(0.,2.*pi,ntheta_out);
rhobound_thetamesh=interpos(thetabndsort,rhobndsort,thetamesh,-0.1,[-1 -1],2.*pi);

nsigma=80; % create a polar rho mesh as fraction of rhobound for each thetamesh
sigma=linspace(0.,1.,nsigma);
for i=1:length(thetamesh)
  rhomesh(1:nsigma,i) = sigma.*rhobound_thetamesh(i);
end

% Compute flux_norm on (rho,theta) mesh by interpolating psi(R,Z) on these Rrho,Rtheta points
thetamesh2D = ones(nsigma,1) * reshape(thetamesh,1,length(thetamesh));
Rrhotheta = Raxis + rhomesh.*cos(thetamesh2D);
Zrhotheta = Zaxis + rhomesh.*sin(thetamesh2D);

tic
[farray_out,varargout]= interpos2Dcartesian(RR,ZZ,farray_in,Rrhotheta,Zrhotheta,tension_default);
toc
% make sure of edge values of normalized flux
farray_out(1,:,1) = 0.;
farray_out(end,:,1) = 1.;

% compute rho polar for each theta corresponding to desired flux surfaces, now that flux_norm known on rho,theta
% use sqrt(flux) since flux not good for inverse interpolation near axis
clear rho_psi
for j=1:ntheta_out
  [rho_psi(:,j)]=interpos(sqrt(farray_out(:,j,1)),rhomesh(:,j),sqrt(flux_norm_out),tension_default,[2 2],[rhomesh(1,j) rhomesh(end,j)]);
  farray_out_fluxnormtheta(:,j,1) = flux_norm_out;
  for k=2:size(farray_out,3)
    farray_out_fluxnormtheta(:,j,k) = interpos(rhomesh(:,j),farray_out(:,j,k),rho_psi(:,j),tension_default, ...
          [2 2],[farray_out(1,j,k) farray_out(end,j,k)]);
  end
end
thetamesh2D = ones(nflux_out,1)*thetamesh;

% rho_psi and farray_out_fluxnormtheta(:,:,k) are now known on (flux_norm_out,thetamesh) mesh
%

% test with plotting flux contours
figure
plot((Raxis + rho_psi.*cos(thetamesh2D))',(Zaxis + rho_psi.*sin(thetamesh2D))')
hold on
contour(RR,ZZ,farray_in(:,:,1)',flux_norm_out,'k--')
break
for k=2:size(farray_out,3)
  figure
  contour(Raxis + rho_psi.*cos(thetamesh2D),Zaxis + rho_psi.*sin(thetamesh2D),farray_out_psitheta(:,:,k),40)
  hold on
  contour(RR,ZZ,farray_in(:,:,k)',40,'k--')
  title([num2str(k)])
end
