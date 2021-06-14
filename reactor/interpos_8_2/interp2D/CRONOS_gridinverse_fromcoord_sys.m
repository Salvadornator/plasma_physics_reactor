load equil;

coord_sys=equil.coord_sys;
psiin1D=coord_sys.grid.dim1; % psi mesh of R(psi,chi), Z(psi,chi)
psiin1D_norm=(psiin1D-psiin1D(1))./(psiin1D(end)-psiin1D(1));
rhopsiin1D_norm=sqrt(psiin1D_norm);
rhopsiin1D_norm(end)=1.;
% compute theta-mesh
R0=coord_sys.position.r(1,1);
Z0=coord_sys.position.z(1,1);
ndim1=length(psiin1D);
% theta at axis no meaning
for i=2:ndim1
  rho_coord(i,:)=sqrt((coord_sys.position.r(i,:)-R0).^2 + (coord_sys.position.z(i,:)-Z0).^2);
  theta_coord(i,:)=atan2(coord_sys.position.z(i,:)-Z0,coord_sys.position.r(i,:)-R0);
end
rho_coord(1,:)=0.;
theta_coord(1,:)=theta_coord(2,:);
% use 0 to 2pi
ii=find(theta_coord<0);
theta_coord(ii) = theta_coord(ii) + 2.*pi;

% $$$ % Work on theta mesh increasing between -pi, +pi and get rho_LCFS
% $$$ rhobound=rho_coord(end,:);
% $$$ thetabound=theta_coord(end,:);
% $$$ [thetabound0,isort]=sort(thetabound);
% $$$ rhobound0=rhobound(isort);
% $$$ thetabound=thetabound0;
% $$$ rhobound=rhobound0;

% psi in center not given
ndim2=length(coord_sys.grid.dim2);
% $$$ theta_coord_eff(2:ndim1,1:ndim2) = theta_coord(2:end,isort);
% $$$ rho_coord_eff(2:ndim1,1:ndim2) = rho_coord(2:end,isort);
% $$$ for j=1:ndim2
% $$$   sigma_rho_coord_eff(:,j) = rho_coord_eff(:,j)./rhobound(j);
% $$$ end
% $$$ sigma_rho_coord_eff(end,:) = 1.;

% psi and theta on which R, Z are required
nrhopsiout=51;
rhopsiout_norm=linspace(rhopsiin1D_norm(1),rhopsiin1D_norm(end),nrhopsiout);
nthetaout=61;
theta_out=linspace(-pi,pi,nthetaout);

% first compute values on each theta_out
tension_theta=-0.1;
for i=2:ndim1
  [thetasorted_i,isort_i]=sort(theta_coord(i,:));
  rho_thetaout(i,:)=interp1(thetasorted_i,rho_coord(i,isort_i),theta_out);
%  rho_thetaout(i,:)=interpos(thetasorted_i,rho_coord(i,isort_i),theta_out,tension_theta,[-1],2.*pi);
end
rho_thetaout(1,:)=0.;
R1=R0+rho_thetaout*cos(theta_out);
Z1=Z0+rho_thetaout*sin(theta_out);
break
% Then compute on psi_out
for j=1:nthetaout
  % relative to psi value
  rho_rhopsiout_thetaout(1:nrhopsiout,j) = interp1(psiin1D_norm,rho_thetaout(:,j),rhopsiout_norm);
  % rho_rhopsiout_thetaout(1:nrhopsiout,j) = interpos(psiin1D_norm,rho_thetaout(:,j),rhopsiout_norm,tension_rho,[2 2],[0. 1.]);
  rho_rhopsiout_thetaout(1:nrhopsiout,j) = rho_rhopsiout_thetaout(1:nrhopsiout,j) .* rhobound_thetaout(j);
  RRout(1:nrhopsiout,j) = R0 + rho_rhopsiout_thetaout(1:nrhopsiout,j) .* cos(theta_out(j));
  ZZout(1:nrhopsiout,j) = Z0 + rho_rhopsiout_thetaout(1:nrhopsiout,j) .* sin(theta_out(j));
end

break

tension_rho= -0.1;
tension_bound = -0.1;
iextrapol=0;
nbc_rho=[1 2];
