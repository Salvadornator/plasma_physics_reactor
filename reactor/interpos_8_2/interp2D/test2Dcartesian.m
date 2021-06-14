set_defaults_matlab

load equil.mat

RR=equil.profiles_2d.grid.dim1;
ZZ=equil.profiles_2d.grid.dim2;
psi_rz=equil.profiles_2d.psi;

Rout=linspace(RR(1),RR(end),101)'*ones(1,111);
Zout=ones(101,1)*linspace(ZZ(1),ZZ(end),111);

farray_in(:,:,1) = equil.profiles_2d.psi;
farray_in(:,:,2) = equil.profiles_2d.br;
farray_in(:,:,3) = equil.profiles_2d.bz;
farray_in(:,:,4) = equil.profiles_2d.bphi;
[farray_out,dfarrdx_out,dfarrdy_out,d2farrdx2_out,d2farrdy2_out,d2farrdxdy_out,varargout]= ...
    interpos2Dcartesian(RR,ZZ,farray_in,Rout,Zout,[],[],[],[],[],[]);

figure
contour(Rout,Zout,farray_out(:,:,1),60);
axis equal

figure
contour(Rout,Zout,farray_out(:,:,2),60);
axis equal

figure
contour(Rout,Zout,farray_out(:,:,3),60);
axis equal

figure
contour(Rout,Zout,farray_out(:,:,4),60);
axis equal

% Compute a rho,theta mesh
Rbnd=equil.eqgeometry.boundary.r;
Zbnd=equil.eqgeometry.boundary.z;

Raxis = equil.global_param.mag_axis.position.r;
Zaxis = equil.global_param.mag_axis.position.z;
psi_axis = equil.global_param.psi_ax;
psi_edge = equil.global_param.psi_bound;

ntheta=201;
theta = linspace(0,2*pi,ntheta);

rhobnd=sqrt((Rbnd-Raxis).^2 + (Zbnd-Zaxis).^2);
thetabnd=atan2(Zbnd-Zaxis,Rbnd-Raxis);
ii=find(thetabnd<0);
thetabnd(ii) = thetabnd(ii) + 2*pi;
[thetasort,isort] = sort(thetabnd);
rhobndsort = rhobnd(isort);

figure
plot(Rbnd,Zbnd)
hold on
plot(Raxis+rhobndsort.*cos(thetasort),Zaxis+rhobndsort.*sin(thetasort),'r*')

figure
plot(thetasort,rhobndsort)

% theta mesh
ntheta=100;
thetamesh=linspace(0.,2.*pi,ntheta+1);
thetamesh = thetamesh(1:end-1); % remove redundant point at 2pi
rhobound=interpos(thetasort,rhobndsort,thetamesh,-0.1,[-1 -1],2.*pi);
nsigma=60;
sigma=linspace(0.,1.,nsigma);
for i=1:length(thetamesh)
  rhomesh(1:nsigma,i) = sigma.*rhobound(i);
end

% Compute psi on (rho,theta) mesh
thetaeff = ones(nsigma,1) * reshape(thetamesh,1,length(thetamesh));
Rrhotheta = Raxis + rhomesh.*cos(thetaeff);
Zrhotheta = Zaxis + rhomesh.*sin(thetaeff);

tic
[farray_out,dfarrdx_out,dfarrdy_out,d2farrdx2_out,d2farrdy2_out,d2farrdxdy_out,varargout]= ...
    interpos2Dcartesian(RR,ZZ,farray_in,Rrhotheta,Zrhotheta,[],[],[],[],[],[]);
toc
farray_out(1,:,1) = psi_axis;
farray_out(end,:,1) = psi_edge;

npsi_out=80;
psi_out=linspace(psi_axis,psi_edge,npsi_out);
dpsisign=sign(psi_edge-psi_axis);

clear rho_psi
for j=1:ntheta
  [rho_psi(:,j)]=interpos(sqrt(dpsisign.*farray_out(:,j,1)),rhomesh(:,j),sqrt(dpsisign.*psi_out),-0.1,[2 2],[rhomesh(1,j) rhomesh(end,j)]);
  for k=2:size(farray_out,3)
    farray_out_psitheta(:,j,k) = interpos(rhomesh(:,j),farray_out(:,j,k),rho_psi(:,j),-0.1);
  end
end
thetamesh2D = ones(npsi_out,1)*thetamesh;
% now we have psi(rho_psi,thetamesh)=farray_out(:,:,1)
figure
plot((Raxis + rho_psi.*cos(thetamesh2D))',(Zaxis + rho_psi.*sin(thetamesh2D))')
hold on
contour(RR,ZZ,farray_in(:,:,1)',psi_out,'k--')

for k=2:size(farray_out,3)
  figure
  contour(Raxis + rho_psi.*cos(thetamesh2D),Zaxis + rho_psi.*sin(thetamesh2D),farray_out_psitheta(:,:,k),40)
  hold on
  contour(RR,ZZ,farray_in(:,:,k)',40,'k--')
  title([num2str(k)])
end
