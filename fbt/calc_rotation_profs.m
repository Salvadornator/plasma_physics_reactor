function out = calc_rotation_profs(g,Ti,d_i,d_e,Rx,Te,ni,nc,Zeff,pion,pc,pathname)

%% Calculating kv-constant for neoclassical poloidal velocity
mi = 1.67e-27;
Zi = 1;
ZI = 6;
e = 1.6e-19;
mI = 12*mi;
e0 = 8.85e-12;
alfa = ZI^2*nc./ni/Zi^2;
beta = (27*mi/4/mI)^2./(15/2 + sqrt(2*alfa));
G = g.ft./g.fc;
v_th = sqrt(e*Ti./mi);
omega_ta = v_th./g.rcentr./g.qpsi;
tau_ii = 12*sqrt(pi^3)*e0^2*mi^2*v_th.^3./(ni*e^4.*coulomb_log_i(ni,Ti,Zeff));
nu_star_ii = 1./((g.amins./g.Rgeo).^1.5.*omega_ta.*tau_ii);

mu00_B_hat = 0.53+Zeff;
K01_B_hat = 0.71+Zeff;
mu11_B_hat = 1.39+13*Zeff/4;

mu00_P_hat = 3.54;
K01_P_hat = 10.63;
mu11_P_hat = 11.52;

Deff = 2.23 + 5.32*Zeff + 2.4*Zeff.^2;
mu00_PS_hat = (3.02+4.25*Zeff)./Deff;
K01_PS_hat = (12.43+20.13*Zeff)./Deff;
mu11_PS_hat = (15.38+26.97*Zeff)./Deff;

mu00_i = G.*mu00_B_hat./(1+2.92*nu_star_ii.*mu00_B_hat./mu00_P_hat)./(1 + mu00_P_hat./(6*omega_ta.*tau_ii.*mu00_PS_hat));
mu11_i = G.*mu11_B_hat./(1+2.92*nu_star_ii.*mu11_B_hat./mu11_P_hat)./(1 + mu11_P_hat./(6*omega_ta.*tau_ii.*mu11_PS_hat));
K01    = G.*K01_B_hat./(1+2.92*nu_star_ii.*K01_B_hat./K01_P_hat)./(1 + K01_P_hat./(6*omega_ta.*tau_ii.*K01_PS_hat));
mu01_i = 5/2*mu00_i - K01;
D = mu00_i.*(mu11_i + sqrt(2) + alfa - alfa.*beta) - mu01_i.^2;
k1 = mu01_i.*(sqrt(2) + alfa - alfa.*beta)./D;
k2 = (mu00_i.*mu11_i - mu01_i.^2)./D;

%% Calculating hydrogen poloidal rotation
dTidpsi = dFdx(linspace(g.simag,g.sibry,length(g.psin)),Ti).';
dpidpsi = dFdx(linspace(g.simag,g.sibry,length(g.psin)),pion).';
dpcdpsi = dFdx(linspace(g.simag,g.sibry,length(g.psin)),pc).';
%out.u_hat = sign(g.current)*k1.*g.fpol./g.B2.*dTidpsi; % Main Ion
out.u_hat = sign(g.current).*g.fpol.*Ti./g.B2*Zi/ZI.*((k1+3/2*k2)/Zi.*dTidpsi./Ti - dpidpsi/ZI./pion + dpcdpsi/ZI./pc); % Carbon

%% Calculating hydrogen toroidal rotation
out.psin = g.psin;
out.rho_tor = g.rho_tor;
R_sep = linspace(g.rmaxis,g.rmaxis+2,10001)';
Z_sep = g.zmaxis*ones(size(R_sep));
F = interp2(g.rg,g.zg,g.psirzn,R_sep,Z_sep,'cubic')'; F(1) = 0;
R_sep = R_sep(F<=1);
R_pt = max(R_sep(F<=1-d_i));
R_sep = max(R_sep);
Lte = (R_sep - R_pt)*(interp1(g.psin,Te,1-d_e,'pchip') + Te(end))/(interp1(g.psin,Te,1-d_e,'pchip') - Te(end))/2;
dc = 0;
if isempty(Rx)
    Rx = g.rcentr;
end
Rxm = (2*Rx - max(g.rbbbs) - min(g.rbbbs))/(max(g.rbbbs) - min(g.rbbbs));
uphi_ped = 1e3*sign(g.current)*0.104*interp1(g.psin,g.fc,1-d_i,'pchip')*(dc/2 - Rxm)*g.q95/(Lte*1e2)*interp1(g.psin,Ti,1-d_i,'pchip')/abs(g.bcentr);
uphi_core = sign(g.current)*(30+0.1*(g.wp-630)/1e3/(abs(g.current)/1e6))*1e3;
uphi = tanhped(g.psin,0,1-d_i/2,d_i,uphi_ped,uphi_core-uphi_ped,1,2);
out.omega_hat = (uphi - out.u_hat.*g.Bt_mp)./g.R_mp;

%% ExB rotation
%dpcdpsi = dFdx(linspace(g.simag,g.sibry,length(g.psin)),pc).';
out.omega_exb = dpcdpsi./nc/1.6e-19/ZI + out.omega_hat;
out.Er = out.omega_exb.*g.R_mp.*g.Bp_mp;

%% Poloidal and toroidal velocities
out.utor = out.omega_hat.*g.R_mp + out.u_hat.*g.Bt_mp;
out.upol = out.u_hat.*g.Bp_mp;

%% Saving profiles
utor = [g.psin out.utor/1e3];
upol = [g.psin out.upol/1e3];
omega_exb = [g.psin out.omega_exb/1e3];
omega_hat = [g.psin out.omega_hat/1e3];
u_hat = [g.psin out.u_hat/1e3];
save(fullfile(pathname,'profile_upol'),'upol','-ascii')
save(fullfile(pathname,'profile_utor'),'utor','-ascii')
save(fullfile(pathname,'profile_omega_exb'),'omega_exb','-ascii')
save(fullfile(pathname,'profile_omega_hat'),'omega_hat','-ascii')
save(fullfile(pathname,'profile_u_hat'),'u_hat','-ascii')

%% Plotting profiles
if 0
    figure(17)
    clf
    plot(g.psin,out.upol_hat/1e3,'r','linewidth',2)
    hold on
    plot(g.psin,out.utor_hat/1e3,'b','linewidth',2)
    xlabel('\psi_N')
    title('Plasma Main Ion Rotation ( krad/s ) - Positive if clockwise')
    legend('v_{\theta}','v_{\phi}','location','northwest')
end

end