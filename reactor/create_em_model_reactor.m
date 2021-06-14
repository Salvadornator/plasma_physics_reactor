function create_em_model_reactor

fprintf('\n       --- Calculating REACTOR static data ---\n')
%% Loading TCABR coils
psirz  = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/greenfunc/green_table_65x65_fbt.mat');
cname  = {'OH' 'E1' 'E2' 'E3' 'E4' 'E5' 'E6' 'E7' 'E8' 'E9' 'E10' 'F1' 'F2' 'F3' 'F4' 'D1' 'D2' 'D3' 'D4'};
Coil   = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_coils');
C = Coil([1 8:end],:);

%% Active-active coil coupling
% Resistance Matrix
Ra = 1.7e-8*Coil(:,5)*2*pi.*Coil(:,1)./(Coil(:,3).*Coil(:,4)./Coil(:,5)); %Surface area divided by the number of turns
Ra = [sum(Ra(1:7)); Ra(8:end)];

%Inductance matrix
bypass = C(1,4);
C(1,4)=0;
L = self_inductance(C(:,1),C(:,3).*C(:,4),C(:,5));
C(1,4)=bypass;
L(1,1) = 0;
for ioh = 1:7
    L(1,1) = L(1,1) + max(interp2(psirz.r,psirz.z,Coil(ioh,5)*psirz.OH.G,linspace(Coil(ioh,1)*0.8,Coil(ioh,1)*1.2,1001),Coil(ioh,2)*ones(1,1001),'cubic'));
end
Maa = diag(L);
for icoil = 1:size(C,1)
    for jcoil = icoil+1:size(C,1)
        eval(['Maa(icoil,jcoil) = mutual_inductance(psirz.r,psirz.z,psirz.' cname{icoil} '.G,C(jcoil,1),C(jcoil,2),C(jcoil,3),C(jcoil,4),C(jcoil,5),10);'])
    end
end
Maa = triu(Maa) + triu(Maa,1)';

%% Vessel-vessel filament coupling
[W,dRv,dZv] = vessel_filaments_reactor;
dA = abs(dRv.*dZv);

% Resistance vector
Rv = 7.4e-7*2*pi*W(:,1)./dA;

%Inductance matrix
Mvv = zeros(size(W,1));
for icoil = 1:size(W,1)
    Mvv(icoil,:) = green_em(W(icoil,1),W(icoil,2),W(:,1),W(:,2),'mutual').';
    Mvv(icoil,icoil) = self_inductance(W(icoil,1),dA(icoil),1);
end

% Reduced number of eigenmodes
[V,D] = eig(-Mvv\diag(Rv));
[~,ind] = sort(diag(D),'descend');
V_e_v = V(:,ind);
for ii = 1:size(V_e_v,2)
    if sum(V_e_v(:,ii))^2 < 1e-24
        V_e_v(:,ii) = V_e_v(:,ii);
    else
        V_e_v(:,ii) = V_e_v(:,ii)/sum(V_e_v(:,ii));
    end
end
T_e_v = inv(V_e_v);

%% Coil-Vessel coupling
Mav = zeros(size(C,1),size(W,1));
Mav(1,:) = (interp2(psirz.r,psirz.z,0.99*psirz.OH.G,W(:,1),W(:,2),'cubic')).';
for icoil = 2:size(C,1)
    eval(['Mav(icoil,:) = (interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.G,W(:,1),W(:,2),''cubic'')).'';']);
end

%% Grid-grid
[R,Z] = meshgrid(psirz.r,psirz.z);
rx = R(:);
zx = Z(:);
Mxx    = zeros(length(rx));
dRMxx  = zeros(length(rx));
dZMxx  = zeros(length(rx));
br_x_x = zeros(length(rx));
bz_x_x = zeros(length(rx));
for ip = 1:length(rx)
    Mxx(ip,:)    = green_em(rx(ip),zx(ip),rx,zx,'mutual').';
    dRMxx(ip,:)  = green_em(rx(ip),zx(ip),rx,zx,'dMdr').';
    dZMxx(ip,:)  = green_em(rx(ip),zx(ip),rx,zx,'dMdz').';
    br_x_x(ip,:) = green_em(rx(ip),zx(ip),rx,zx,'br').';
    bz_x_x(ip,:) = green_em(rx(ip),zx(ip),rx,zx,'bz').';
end

%% Grid-coils mutuals and fields
Max       = zeros(size(C,1),length(rx));
dRMax     = zeros(size(C,1),length(rx));
dZMax     = zeros(size(C,1),length(rx));
br_a_x    = zeros(size(C,1),length(rx));
bz_a_x    = zeros(size(C,1),length(rx));
dbrdr_a_x = zeros(size(C,1),length(rx));
dbzdr_a_x = zeros(size(C,1),length(rx));
dbrdz_a_x = zeros(size(C,1),length(rx));
dbzdz_a_x = zeros(size(C,1),length(rx));
for icoil = 1:size(C,1)
    eval(['Max(icoil,:)       = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.G,    rx,zx,''cubic'').'';']);
    eval(['dRMax(icoil,:)     = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dMdR, rx,zx,''cubic'').'';']);
    eval(['dZMax(icoil,:)     = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dMdZ, rx,zx,''cubic'').'';']);
    eval(['br_a_x(icoil,:)    = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.BR,   rx,zx,''cubic'').'';']);
    eval(['bz_a_x(icoil,:)    = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.BZ,   rx,zx,''cubic'').'';']);
    eval(['dbrdr_a_x(icoil,:) = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dBRdR,rx,zx,''cubic'').'';']);
    eval(['dbzdr_a_x(icoil,:) = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dBZdR,rx,zx,''cubic'').'';']);
    eval(['dbrdz_a_x(icoil,:) = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dBRdZ,rx,zx,''cubic'').'';']);
    eval(['dbzdz_a_x(icoil,:) = interp2(psirz.r,psirz.z,psirz.' cname{icoil} '.dBZdZ,rx,zx,''cubic'').'';']);
end

%% Grid-vessel
Mvx   = zeros(size(W,1),length(rx));
dRMvx = zeros(size(W,1),length(rx));
dZMvx = zeros(size(W,1),length(rx));
for ivv = 1:size(W,1)
    Mvx(ivv,:)   = green_em(W(ivv,1),W(ivv,2),rx,zx,'mutual').';
    dRMvx(ivv,:) = green_em(W(ivv,1),W(ivv,2),rx,zx,'dMdr').';
    dZMvx(ivv,:) = green_em(W(ivv,1),W(ivv,2),rx,zx,'dMdz').';
end

%% Output
G.rx = reshape(rx,65,65); G.rx = G.rx(1,:)';
G.zx = reshape(zx,65,65); G.zx = G.zx(:,1);
G.rxx = rx;
G.zxx = zx;
G.Mxx = Mxx;
G.dRMxx = dRMxx;
G.dZMxx = dZMxx;
G.bz_x_x = bz_x_x;
G.br_x_x = br_x_x;
G.Max = Max;
G.dRMax = dRMax;
G.dZMax = dZMax;
G.bz_a_x = bz_a_x;
G.br_a_x = br_a_x;
G.dbzdr_a_x = dbzdr_a_x;
G.dbrdr_a_x = dbrdr_a_x;
G.dbzdz_a_x = dbzdz_a_x;
G.dbrdz_a_x = dbrdz_a_x;
G.Mxa = Max.';
G.dRMxa = dRMax.';
G.dZMxa = dZMax.';
G.bz_x_a = bz_a_x.';
G.br_x_a = br_a_x.';
G.dbzdr_x_a = dbzdr_a_x.';
G.dbrdr_x_a = dbrdr_a_x.';
G.dbzdz_x_a = dbzdz_a_x.';
G.dbrdz_x_a = dbrdz_a_x.';
G.Maa  = Maa;
G.Mvv  = Mvv;
G.Mvx  = Mvx;
G.dRMvx = dRMvx;
G.dZMvx = dZMvx;
G.Mxv  = Mvx.';
G.dRMxv = dRMvx.';
G.dZMxv = dZMvx.';
G.Mav  = Mav;
G.Mva  = Mav.';
G.Ra   = Ra;
G.Rv   = Rv;
G.V_e_v = V_e_v;
G.T_e_v = T_e_v;
G.V    = V;
G.labels_a = cname;

%% Saving TCABR STATIC data
fprintf('       --- Saving REACTOR static data ---\n')
save('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_static','-struct','G')
fprintf('       --- STATIC data saved: /home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_static\n')

end