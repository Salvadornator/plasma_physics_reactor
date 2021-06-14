function params = j2psi(params,iplot)

%% Calculating Green's functions between plasma filaments and grid points
fprintf('       --- Calculating poloidal flux due to the plasma\n')
Jp = params.jcur.jplasma;
r = params.jcur.rg;
z = params.jcur.zg;
dR = mean(diff(r));
dZ = mean(diff(z));

[Rg,Zg] = meshgrid(r,z);
Rgg = Rg(:);
Zgg = Zg(:);
Jp = Jp(:);
isin = find(inpolygon(Rgg,Zgg,params.Rbb,params.Zbb));
G = zeros(length(z),length(r)); BR = G; BZ = G; 
if strcmp(params.plasmaid,'snowflake')
    dBRdR = G; dBRdZ = G; dBZdR = G; dBZdZ = G;
end

for iplasma = 1:length(isin)
        G  = G + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'psi');
        BR = BR + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'br');
        BZ = BZ + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'bz');
        if strcmp(params.plasmaid,['SF','snowflake'])
            dBRdR = dBRdR + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'dbrdr');
        	dBRdZ = dBRdZ + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'dbrdz');
        	dBZdR = dBZdR + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'dbzdr');
        	dBZdZ = dBZdZ + Jp(isin(iplasma))*dR*dZ*green_em(Rgg(isin(iplasma)),Zgg(isin(iplasma)),Rg,Zg,'dbzdz');
        end
end
%keyboard
%G = params.Mxx.*J;

%% Interpolating into original meshgrid
[Rg,Zg] = meshgrid(params.rmesh,params.zmesh);
params.psip = interp2(r,z,-G/2/pi,Rg,Zg,'cubic');
params.brp = interp2(r,z,BR,Rg,Zg,'cubic');
params.bzp = interp2(r,z,BZ,Rg,Zg,'cubic');
if strcmp(params.plasmaid,['SF','snowflake'])
    params.greens.dbrdrp = interp2(r,z,dBRdR,Rg,Zg,'cubic');
    params.greens.dbrdzp = interp2(r,z,dBRdZ,Rg,Zg,'cubic');
    params.greens.dbzdrp = interp2(r,z,dBZdR,Rg,Zg,'cubic');
    params.greens.dbzdzp = interp2(r,z,dBZdZ,Rg,Zg,'cubic');
end

%% Finding magnetic axis
[~,ind] = min(G(:));
Rmaxis = Rg(:); Rmaxis = Rmaxis(ind);
Zmaxis = Zg(:); Zmaxis = Zmaxis(ind);

%% Output
params.jcur.G = -G/2/pi;
params.jcur.BR = BR;
params.jcur.BZ = BZ;
if strcmp(params.plasmaid,['SF','snowflake'])
    params.jcur.dBRdR = dBRdR;
    params.jcur.dBZdR = dBZdR;
    params.jcur.dBRdZ = dBRdZ;
    params.jcur.dBZdZ = dBZdZ;
end
params.jcur.rmaxis = Rmaxis;
params.jcur.zmaxis = Zmaxis;

%% Plotting J distribution
if iplot
    figure(1)
    grid on
    clf
    contourf(params.rmesh,params.zmesh,params.psip,21)
    hold on
    plot(params.Rbb,params.Zbb,'r','linewidth',3)
    plot(params.Rmaxis,params.Zmaxis,'w+','linewidth',3,'markersize',10)
    plot(params.Rxp,params.Zxp,'wx','linewidth',3,'markersize',10)
    plot_tcabr('fig',1)
    hold off
    title('Plasma Poloidal Flux Distribution')
    %caxis([0 ceil(max(params.jcur.jplasma(:)/1e6))-1])
    colorbar
    axis([0.35 0.9 -0.3001 0.3001])
    drawnow
    %%
end

end