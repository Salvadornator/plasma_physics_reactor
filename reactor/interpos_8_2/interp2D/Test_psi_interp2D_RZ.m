
set_defaults_matlab

load for_griddata.mat
load rho_griddata.mat
load psi_griddata

ndim1=size(reqk,1);
ndim2=size(reqk,2);
theta=linspace(0,2.*pi,ndim2+1);
for i=1:ndim1
  rho_eqk(i,:)=sqrt((reqk(i,:)-reqk(1,1)).^2 + (zeqk(i,:)-zeqk(1,1)).^2);
  theta_eqk(i,:)=atan2(zeqk(i,:)-zeqk(1,1),reqk(i,:)-reqk(1,1));
end
% It seems ndim2 is a theta mesh from 0 to 2pi
% extract theta_mesh
ij=find(theta_eqk<-1e-6);
thet2=theta_eqk;
thet2(ij)=theta_eqk(ij)+2*pi;
for j=1:ndim2
  theta_mesh(j)=mean(thet2(2:end,j));
end

tension1D_R = -0.1;
tension1D_Z = -0.1;

ioptos=13;

% we have psi(rhoRZ,theta) in psiRZ and we want to obtain it on RRgrid, ZZgrid
% R0,Z0 of rho,theta mesh is reqk(1,1), zeqk(1,1)
% Do the interpolation/extrapolation on rho,theta "space" to use periodic conditions and extrapolate along theta=cst
%
% Find rho, theta values at RRgrid,ZZgrid
for i=1:size(RRgrid,1)
  rhoRRgrid(i,:) = sqrt((RRgrid(i,:)-reqk(1,1)).^2 + (ZZgrid(i,:)-zeqk(1,1)).^2);
  thetaRRgrid(i,:) = atan2(ZZgrid(i,:)-zeqk(1,1),RRgrid(i,:)-reqk(1,1));
end

% We need a 2D interpolation for each "line" in (rho,theta) plane
tension_rho= -0.1;
tension_theta=-0.1;
iextrapol=1;
nbc_rho=[2 2];
rhoRZ=reshape(rhoRZ,length(rhoRZ),1);
theta=reshape(theta,1,length(theta));
psi_in=reshape(psiRZ,length(psiRZ),1) * ones(1,length(theta));
break
[psi_out,drho,dtheta,d2rho,d2theta,d2rhotheta]=interpos2Dpolar(rhoRZ,theta,psi_in,rhoRRgrid,thetaRRgrid,tension_rho,tension_theta,nbc_rho,iextrapol);

