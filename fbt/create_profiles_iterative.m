function struct = create_profiles_iterative(ne0,Te0,Ti0,nesep,Tesep,Teped,Tiped,neped,d_e,d_i,Zeff,fz,pathname)

%% Calculating profiles
struct.psin = linspace(0,1.0,129)';
struct.Zeff = Zeff*ones(size(struct.psin));
if fz < (Zeff-1)/650
    struct.fz = fz*ones(size(struct.psin));
else
    disp(['   --- The value of fz was changed from ' num2str(fz*100) '% to ' num2str((Zeff-1)/650*100) '%'])
    fz = (Zeff-1)/650;
    struct.fz = fz*ones(size(struct.psin));
end
struct.ne   = 1e20*tanhped(struct.psin,0,1 - d_e/2,d_e,neped,ne0-neped,2,3);
struct.ne   = struct.ne - struct.ne(end) + nesep*1e20;
struct.Te   = 1e3*tanhped(struct.psin,0,1 - d_e/2,d_e,Teped,Te0-Teped,2,3);
struct.Te   = struct.Te - struct.Te(end) + Tesep*1e3;
struct.Ti   = 1e3*tanhped(struct.psin,0,1 - d_i/2,d_i,Tiped,Ti0-Tiped,2,3);
struct.Ti   = struct.Ti - struct.Ti(end) + Tesep*1e3;
struct.ni   = struct.ne.*(520*struct.fz+6-struct.Zeff)/5;
struct.nc   = (struct.Zeff-1-650*struct.fz).*struct.ne/30;
struct.nFe  = struct.ne.*struct.fz;
struct.pe   = 16.02e3*struct.ne.*struct.Te/1e20/1e3;
struct.pion = 16.02e3*struct.ni.*struct.Ti/1e20/1e3;
struct.pc   = 16.02e3*struct.nc.*struct.Ti/1e20/1e3;
struct.pFe  = 16.02e3*struct.nFe.*struct.Ti/1e20/1e3;
struct.ptot = struct.pe + struct.pion + struct.pc + struct.pFe;

%% Saving profiles
if 1
    ptot = [struct.psin struct.ptot/1e3];
    ti = [struct.psin struct.Ti/1e3];
    te = [struct.psin struct.Te/1e3];
    ne = [struct.psin struct.ne/1e20];
    ni = [struct.psin struct.ni/1e20];
    nc = [struct.psin struct.nc/1e20];
    nFe = [struct.psin struct.nFe/1e20];
    zeff = [struct.psin struct.Zeff];
    fz = [struct.psin struct.fz];
    save(fullfile(pathname,'profile_p'),'ptot','-ascii')
    save(fullfile(pathname,'profile_ti'),'ti','-ascii')
    save(fullfile(pathname,'profile_te'),'te','-ascii')
    save(fullfile(pathname,'profile_ne'),'ne','-ascii')
    save(fullfile(pathname,'profile_ni'),'ni','-ascii')
    save(fullfile(pathname,'profile_nc'),'nc','-ascii')
    save(fullfile(pathname,'profile_nFe'),'nFe','-ascii')
    save(fullfile(pathname,'profile_zeff'),'zeff','-ascii')
    save(fullfile(pathname,'profile_fz'),'fz','-ascii')
    fprintf('\n')
    disp(['   --- Kinetic profiles SAVED at ' pathname])
else
    fprintf('\n')
    disp('   --- Kinetic profiles NOT SAVED')
end

%% Plotting profiles
figure(106)
clf
hs(1) = subplot(1,3,1);
plot(struct.psin,struct.pe/1e3,'r','linewidth',2)
hold on
plot(struct.psin,struct.pion/1e3,'b','linewidth',2)
plot(struct.psin,struct.pc/1e3,'color',[1 1 1]*0.8,'linewidth',2)
plot(struct.psin,struct.pFe/1e3,'color',[1 1 1]*0.5,'linewidth',2)
plot(struct.psin,struct.ptot/1e3,'k','linewidth',2)
xlabel('\psi_N')
title('Plasma Pressure ( kPa )')
legend('p_e','p_i','p_C','p_{Fe}','p_{total}')

hs(2) = subplot(1,3,2);
plot(struct.psin,struct.Te,'r','linewidth',2)
hold on
plot(struct.psin,struct.Ti,'b','linewidth',2)
xlabel('\psi_N')
title('Temperature ( eV )')
legend('T_e','T_i')

hs(3) = subplot(1,3,3);
plot(struct.psin,struct.ne/1e19,'r','linewidth',2)
hold on
plot(struct.psin,struct.ni/1e19,'b','linewidth',2)
plot(struct.psin,6*struct.nc/1e19,'k','linewidth',2)
plot(struct.psin,26*struct.nFe/1e19,'color',[1 1 1]*0.5,'linewidth',2)
xlabel('\psi_N')
title('Number Density ( 1\times10^{19} m^{-3} )')
legend({'n_e','n_i','6 n_C','26 n_{Fe}'},'location','northeast')

figure(400)
plot(struct.psin,struct.ptot/1e3,'k','linewidth',2)
xlabel('\psi_N')
title('Plasma Total Pressure (p_{tot}) ( kPa )')
xlim([.85,1])
annotation('textbox',[.3 .25 .33 .35],'String','p_{ped}','EdgeColor','none')
annotation('textbox',[.5 .1 .6 .15],'String','\Delta_{\psi}','EdgeColor','none')
linkaxes(hs,'x')

end