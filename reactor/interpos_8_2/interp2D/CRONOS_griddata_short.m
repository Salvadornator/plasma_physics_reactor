load for_griddata.mat

% uses reqk, zeqk, psiold as input

ndim1=size(reqk,1);
ndim2=size(reqk,2);
for i=1:ndim1
  rho_eqk(i,:)=sqrt((reqk(i,:)-reqk(1,1)).^2 + (zeqk(i,:)-zeqk(1,1)).^2);
  theta_eqk(i,:)=atan2(zeqk(i,:)-zeqk(1,1),reqk(i,:)-reqk(1,1));
end

% Work on theta mesh increasing between -pi, +pi and get rho_LCFS
rhobound=rho_eqk(end,:);
thetabound=theta_eqk(end,:);
[thetabound0,isort]=sort(thetabound);
rhobound0=rhobound(isort);
% $$$ thetabound=linspace(-pi,pi,length(thetabound)+1);
% $$$ thetabound=thetabound(1:end-1);
thetabound=thetabound0;
rhobound=rhobound0;

% psi in center not given
theta_eqk_eff = theta_eqk(2:end,isort);
rho_eqk_eff = rho_eqk(2:end,isort);
psi_in_eff = double(psiold(:,isort));

tension_rho= -0.1;
tension_theta=-0.1;
tension_bound = -0.1;
iextrapol=0;
nbc_rho=[1 2];

% sigma mesh: Use the fact that the input psiold is given at fixed theta values to interpolate along rho on fixed sigma values
sigma=linspace(0,1,100);
clear rho_sigma_eqk psi_in_sigma_theta
for itheta_eqk=1:length(thetabound)
  rho_sigma_eqk(:,itheta_eqk)=sigma.*rho_eqk_eff(end,itheta_eqk);
  psi_in_sigma_theta(:,itheta_eqk)=interpos(rho_eqk_eff(:,itheta_eqk),psi_in_eff(:,itheta_eqk),rho_sigma_eqk(:,itheta_eqk),tension_rho,[1 2],[0. psi_in_eff(end,itheta_eqk)]);
end

testRZgrids=0;
switch testRZgrids
 case 0
  % uses input RRgrid, ZZgrid
 case 1
  Rgrid=linspace(1.5,4.5,100)';
  Zgrid=linspace(-2,2,200);
  RRgrid=Rgrid * ones(size(Zgrid));
  ZZgrid=ones(size(Rgrid)) * Zgrid;
 case 2
  % one chord
  RZstart=[4,-1.5];
  RZend=[2,1.5];
  nbpoints=120;
  RRgrid=linspace(RZstart(1),RZend(1),nbpoints)';
  ZZgrid=linspace(RZstart(2),RZend(2),nbpoints)';
end

clear rhoRRgrid thetaRRgrid rhobound_thetaRRgrid sigma_rhoRRgrid psi_out
for i=1:size(RRgrid,1)
  rhoRRgrid(i,:) = sqrt((RRgrid(i,:)-reqk(1,1)).^2 + (ZZgrid(i,:)-zeqk(1,1)).^2);
  thetaRRgrid(i,:) = atan2(ZZgrid(i,:)-zeqk(1,1),RRgrid(i,:)-reqk(1,1));
  rhobound_thetaRRgrid(i,:) = interpos(thetabound,rhobound,thetaRRgrid(i,:),tension_bound,[-1],2*pi);
  sigma_rhoRRgrid(i,:) = rhoRRgrid(i,:)./rhobound_thetaRRgrid(i,:);
end
[psi_out,drho,dtheta,d2rho,d2theta,d2rhotheta]=interpos2Dpolar(sigma,thetabound,psi_in_sigma_theta,sigma_rhoRRgrid,thetaRRgrid,tension_rho,tension_theta,nbc_rho,iextrapol);

psiaxis = interpos(reqk(1:end,1),[psiold(1,1); psiold(1:end,1)],reqk(1,1),tension_rho,[1 2],[0. psiold(end,1)],[100 ones(1,size(psiold,1))]);
psiedge = psiold(end,1);
if testRZgrids~=2
  figure
  % psivals=linspace(psiaxis,2.*psiedge-psiaxis,100);
  contour(RRgrid,ZZgrid,psi_out,100)
  hold on
  plot(reqk(1,1)+rhobound.*cos(thetabound),zeqk(1,1)+rhobound.*sin(thetabound),'k-')
  axis equal
else
  figure
  contour(reqk(2:end,:),zeqk(2:end,:),psiold,100);
  hold on
  plot([RZstart(1) RZend(1)],[RZstart(2) RZend(2)],'k')
  plot(reqk(1,1)+rhobound.*cos(thetabound),zeqk(1,1)+rhobound.*sin(thetabound),'k-')
  axis equal
  
  figure
  ij=find(psi_out>=psiedge);
  subplot(2,1,1)
  plot(RRgrid,psi_out,'*-')
  hold on
  plot(RRgrid(ij),psi_out(ij),'r-')
  xlabel('R')
  title('in red part inside plasma')
  subplot(2,1,2)
  plot(ZZgrid,psi_out,'*-')
  hold on
  plot(ZZgrid(ij),psi_out(ij),'r-')
  xlabel('Z')
end
