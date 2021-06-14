function out = flux_average(out)

%% Calculating magnetic shear
% Calculating magnetic fields
Muo = pi*4e-7;
dr = mean(diff(out.rg));
dz = mean(diff(out.zg));
[R,Z] = meshgrid(out.rg,out.zg);
out.R = R;
out.Z = Z;
ind = inpolygon(R,Z,out.rbbbs,out.zbbbs);
[dpsidr,dpsidz] = gradient(out.psirz,dr,dz);
out.dpsirzdr = dpsidr;
out.dpsirzdz = dpsidz;
pprime = reshape(interp1(out.psin,out.pprime,out.psirzn(:)),out.nh,out.nh); pprime(~ind) = 0;
out.pprimerz = pprime;
ffprime = reshape(interp1(out.psin,out.ffprim,out.psirzn(:)),out.nh,out.nh); ffprime(~ind) = 0;
out.ffprimrz = ffprime;
fpol = reshape(interp1(out.psin,out.fpol,out.psirzn(:)),out.nh,out.nh); fpol(~ind) = out.fpol(end);
out.fpolrz = fpol;
[dfpoldr,dfpoldz] = gradient(fpol,dr,dz);
R_us = linspace(out.rmaxis,0.815,10001)';
Z_us = out.zmaxis*ones(size(R_us));
F = interp2(out.rg,out.zg,out.psirzn,R_us,Z_us,'cubic')'; F(1) = 0;
R_us = R_us(F<=1);
Z_us = Z_us(F<=1);
psi_mp = interp2(out.rg,out.zg,out.psirzn,R_us,Z_us,'cubic')';
out.R_mp = interp1(psi_mp,R_us,out.psin,'pchip');
out.r_mp = out.R_mp - out.rmaxis;
out.R_ped95 = max(R_us(F<=0.95));
out.fce = 1.6e-19/9.11e-31*abs(out.fpol(end))./out.R_mp/2/pi;
out.fci = 1.6e-19/1.67e-27*abs(out.fpol(end))./out.R_mp/2/pi;
out.jr = -1/Muo./R.*dfpoldz; out.jr(~ind) = 0;
out.jz =  1/Muo./R.*dfpoldr; out.jz(~ind) = 0;
out.jphi = -sign(out.current)*(R.*pprime + 1/Muo./R.*ffprime); out.jphi(~ind) = 0;
BR =  sign(out.current)./R.*dpsidz;
BZ = -sign(out.current)./R.*dpsidr;
BPOL = sqrt(BR.^2 + BZ.^2);
out.br = BR; out.bz = BZ; out.bphi = fpol./R;
[dbrdr,dbrdz] = gradient(out.br,mean(diff(out.rg)),mean(diff(out.zg)));
[dbzdr,dbzdz] = gradient(out.bz,mean(diff(out.rg)),mean(diff(out.zg)));
out.dbrdr = dbrdr; out.dbrdz = dbrdz;
out.dbzdr = dbzdr; out.dbzdz = dbzdz;
out.Bp_mp = interp2(out.rg,out.zg,sqrt(out.br.^2 + out.bz.^2),out.R_mp,out.zmaxis*ones(size(out.R_mp)),'cubic');
out.Bt_mp = interp2(out.rg,out.zg,out.bphi,out.R_mp,out.zmaxis*ones(size(out.R_mp)),'cubic');
out.Bp_ped95 = interp2(out.rg,out.zg,sqrt(out.br.^2 + out.bz.^2),out.R_ped95,out.zmaxis,'cubic');

% Calculating the Jacobian
out.thetarz = -atan2(out.Z - out.zmaxis,out.R - out.rmaxis);
out.thetarz(out.thetarz<0) = out.thetarz(out.thetarz<0) + 2*pi;
[out.dsinthetadr,out.dsinthetadz] = gradient(sin(out.thetarz),dr,dz);
out.dthetadr = out.dsinthetadr./cos(out.thetarz);
out.dthetadz = out.dsinthetadz./cos(out.thetarz);
out.jacobian = -(out.dpsirzdr.*out.dthetadz - out.dpsirzdz.*out.dthetadr)./out.R;
%thetarz = -atan2(out.Z - out.zmaxis,out.R - out.rmaxis);
%[dthetadr,dthetadz] = gradient(thetarz,dr,dz);
%jacobian2 = -(out.dpsirzdr.*dthetadz - out.dpsirzdz.*dthetadr)./out.R;
%out.jacobian(abs(out.jacobian) >= 4) = jacobian2(abs(out.jacobian) >= 4);

% Getting contour lines
psin = [out.psin(3:end-1); 0.999];
C = contourc(out.rg,out.zg,out.psirzn,psin);
C_out = get_contour_levels(C,out.rmaxis,out.zmaxis);
if isempty(C_out)
    out = [];
    return
end

% Output
out.kappa   = zeros(1,length(C_out.n_points));
out.deltau  = zeros(1,length(C_out.n_points));
out.deltal  = zeros(1,length(C_out.n_points));
out.delta  = zeros(1,length(C_out.n_points));
out.qpsi    = zeros(1,length(C_out.n_points));
out.jpar    = zeros(1,length(C_out.n_points));
out.jtor    = zeros(1,length(C_out.n_points));
out.jpol    = zeros(1,length(C_out.n_points));
out.jefit   = zeros(1,length(C_out.n_points));
out.jelite  = zeros(1,length(C_out.n_points));
out.jeliten = zeros(1,length(C_out.n_points));
out.B       = zeros(1,length(C_out.n_points));
out.Bdl     = zeros(1,length(C_out.n_points));
out.B_1     = zeros(1,length(C_out.n_points));
out.B_2     = zeros(1,length(C_out.n_points));
out.B2      = zeros(1,length(C_out.n_points));
out.Bp2     = zeros(1,length(C_out.n_points));
out.Bphi_2  = zeros(1,length(C_out.n_points));
out.Bp      = zeros(1,length(C_out.n_points));
out.Bp_1    = zeros(1,length(C_out.n_points));
out.Ba      = zeros(1,length(C_out.n_points));
out.Rgeo    = zeros(1,length(C_out.n_points));
out.R02_R2  = zeros(1,length(C_out.n_points));
out.R2_1    = zeros(1,length(C_out.n_points));
out.R0_R    = zeros(1,length(C_out.n_points));
out.Vprime  = zeros(1,length(C_out.n_points));
out.A       = zeros(1,length(C_out.n_points));
out.amins   = zeros(1,length(C_out.n_points));
out.fc      = zeros(1,length(C_out.n_points));
out.ft      = zeros(1,length(C_out.n_points));
out.Lpol    = zeros(1,length(C_out.n_points));
out.Ba      = zeros(1,length(C_out.n_points));
out.h       = zeros(1,length(C_out.n_points));
out.h2      = zeros(1,length(C_out.n_points));
out.hh      = zeros(1,length(C_out.n_points));
out.Sperp   = zeros(1,length(C_out.n_points));

for ii = 1:length(C_out.n_points)
    out.Rgeo(ii)  = (min(C_out.Rb{ii}) + max(C_out.Rb{ii}))/2;
    out.amins(ii) = (max(C_out.Rb{ii}) - min(C_out.Rb{ii}))/2;
    dl = [0 sqrt((C_out.Rb{ii}(2:end) - C_out.Rb{ii}(1:end-1)).^2 + (C_out.Zb{ii}(2:end) - C_out.Zb{ii}(1:end-1)).^2)];
    l  = cumsum(dl);
    out.Lpol(ii) = l(end);
    out.Ba(ii) = 4e-7*pi*out.current/out.Lpol(ii);
    R_zmax = C_out.Rb{ii}; R_zmax = mean(R_zmax(max(C_out.Zb{ii}) == C_out.Zb{ii}));
    R_zmin = C_out.Rb{ii}; R_zmin = mean(R_zmin(min(C_out.Zb{ii}) == C_out.Zb{ii}));
    out.kappa(ii)  = (max(C_out.Zb{ii}) - min(C_out.Zb{ii}))/(max(C_out.Rb{ii}) - min(C_out.Rb{ii}));
    out.deltau(ii) = (max(C_out.Rb{ii}) + min(C_out.Rb{ii}) - 2*R_zmax)/(max(C_out.Rb{ii}) - min(C_out.Rb{ii}));
    out.deltal(ii) = (max(C_out.Rb{ii}) + min(C_out.Rb{ii}) - 2*R_zmin)/(max(C_out.Rb{ii}) - min(C_out.Rb{ii}));

    % Interpolating profiles
    bpol           = interp2(out.rg,out.zg,BPOL,C_out.Rb{ii},C_out.Zb{ii},'cubic');
    bphi           = interp1(linspace(0,1,out.nh),out.fpol,C_out.value(ii),'pchip')./C_out.Rb{ii};
    pprime         = interp1(linspace(0,1,out.nh),out.pprime,C_out.value(ii),'pchip');
    ffprime        = interp1(linspace(0,1,out.nh),out.ffprim,C_out.value(ii),'pchip');
    fpol           = interp1(linspace(0,1,out.nh),out.fpol,C_out.value(ii),'pchip');
    h              = sqrt(bpol.^2 + bphi.^2)./max(sqrt(bpol.^2 + bphi.^2));

    % Averages
    Norm           = trapz(l,1./bpol);
    out.Bdl(ii)    = trapz(l,bpol);
    out.Bp(ii)     = trapz(l,bpol./bpol)/Norm;
    out.Bp2(ii)     = trapz(l,bpol.^2./bpol)/Norm;
    out.Bp_1(ii)   = trapz(l,1./bpol./bpol)/Norm;
    out.Vprime(ii) = 2*pi*Norm;
    out.B(ii)      = trapz(l,sqrt(bpol.^2 + bphi.^2)./bpol)/Norm;
    out.h(ii)      = trapz(l,h./bpol)/Norm;
    out.h2(ii)     = trapz(l,h.^2./bpol)/Norm;
    out.hh(ii)     = trapz(l,(1-sqrt(1-h).*(1+0.5*h))./h.^2./bpol)/Norm;
    out.B_1(ii)    = trapz(l,1./sqrt(bpol.^2 + bphi.^2)./bpol)/Norm;
    out.B_2(ii)    = trapz(l,1./(bpol.^2 + bphi.^2)./bpol)/Norm;
    out.B2(ii)     = trapz(l,(bpol.^2 + bphi.^2)./bpol)/Norm;
    out.Bp2(ii)    = trapz(l,bpol.^2./bpol)/Norm;
    out.Bphi_2(ii) = trapz(l,1./bphi.^2./bpol)/Norm;
    out.R02_R2(ii) = trapz(l,(out.rcentr./C_out.Rb{ii}).^2./bpol)/Norm;
    out.R2_1(ii)   = trapz(l,(1./C_out.Rb{ii}.^2)./bpol)/Norm;
    out.R0_R(ii)   = trapz(l,(out.rcentr./C_out.Rb{ii})./bpol)/Norm;
    out.A(ii)      = polyarea(C_out.Rb{ii},C_out.Zb{ii});
    out.Sperp(ii)  = trapz(l,2*pi*C_out.Rb{ii});
    out.qpsi(ii)   = sign(out.bcentr)*fpol*out.Vprime(ii)/(2*pi)^2*out.R2_1(ii);
    
    % Current densities
    out.jpar(ii)   = -sign(out.current)*(out.rcentr*pprime + ffprime/out.rcentr/4e-7/pi*out.B2(ii)/(fpol/out.rcentr).^2);
    out.jtor(ii)   = -sign(out.current)*(out.rcentr*pprime + ffprime/out.rcentr/4e-7/pi*out.R02_R2(ii));
    out.jpol(ii)   = out.jpar(ii) - out.jtor(ii);
    out.jefit(ii)  = out.jtor(ii)/out.R0_R(ii);
end

% Including psin = 0 and 1 in psin vector
psin   = out.psin;
out.qpsi  = interp1(C_out.value,out.qpsi,psin,'pchip');
out.jpar   = interp1(C_out.value,out.jpar,psin,'pchip');
out.fc     = interp1(C_out.value,out.fc,psin,'pchip');
out.ft     = interp1(C_out.value,out.ft,psin,'pchip');
out.jtor   = interp1(C_out.value,out.jtor,psin,'pchip');
out.jpol   = interp1(C_out.value,out.jpol,psin,'pchip');
out.jefit  = interp1(C_out.value,out.jefit,psin,'pchip');
out.B2     = interp1(C_out.value,out.B2,psin,'pchip'); %out.B2 = abs(out.B2);
out.Bp2    = interp1(C_out.value,out.Bp2,psin,'pchip'); %out.Bp2 = abs(out.Bp2);
out.Bp     = interp1(C_out.value,out.Bp,psin,'pchip');
out.Bp_1   = interp1(C_out.value,out.Bp_1,psin,'pchip');
out.Bphi_2 = interp1(C_out.value,out.Bphi_2,psin,'pchip');
out.Bdl    = interp1(C_out.value,out.Bdl,psin,'pchip');
out.R02_R2 = interp1(C_out.value,out.R02_R2,psin,'pchip'); %out.R02_R2 = abs(out.R02_R2);
out.R0_R   = interp1(C_out.value,out.R0_R,psin,'pchip');
out.R2_1   = interp1(C_out.value,out.R2_1,psin,'pchip'); %out.R2_1 = abs(out.R2_1);
out.Rgeo   = interp1(C_out.value,out.Rgeo,psin,'pchip');
out.B      = interp1(C_out.value,out.B,psin,'pchip');
out.B_1    = interp1(C_out.value,out.B_1,psin,'pchip');
out.B_2    = interp1(C_out.value,out.B_2,psin,'pchip');
out.h      = interp1(C_out.value,out.h,psin,'pchip');
out.h2     = interp1(C_out.value,out.h2,psin,'pchip');
out.hh     = interp1(C_out.value,out.hh,psin,'pchip');
out.kappa  = interp1(C_out.value,out.kappa,psin,'pchip');
out.deltal = interp1(C_out.value,out.deltal,psin,'pchip');
out.deltau = interp1(C_out.value,out.deltau,psin,'pchip');
out.A      = interp1(C_out.value,out.A,psin,'pchip');
out.Sperp  = interp1(C_out.value,out.Sperp,psin,'pchip');
out.amins  = interp1(C_out.value,out.amins,psin,'pchip');
out.Vprime = interp1(C_out.value,out.Vprime,psin,'pchip');
out.Ba     = interp1(C_out.value,out.Ba,psin,'pchip');
out.Lpol   = interp1(C_out.value,out.Lpol,psin,'pchip');
out.delta  = (out.deltal + out.deltau)/2;
out.Vp     = out.Vprime/2/pi;
out.V      = cumtrapz(out.psin,out.Vprime)*(out.sibry - out.simag);
out.alpha  = 2e-7/pi*out.Vprime.*abs(out.pprime).*sqrt(out.V/2/pi/pi/out.rcentr);
out.alpha_max = max(out.alpha(out.psin>0.9));
out.qprime = dFdx(out.psin,out.qpsi)/(out.sibry - out.simag);
out.s      = 2*out.V./out.Vprime.*out.qprime./out.qpsi;
out.phi_tor = 2*pi*cumtrapz(linspace(out.simag,out.sibry,length(out.psin)),out.qpsi);
out.rho_tor = sqrt(out.phi_tor/out.phi_tor(end));
%out.j_ps   = -sign(out.current)*out.rcentr*out.pprime.*(1 - out.B2.*out.Bphi_2);
out.wp     = 3/2*trapz(out.V,out.pres);
out.betap  = trapz(out.V,out.pres)/out.V(end)*2*4e-7*pi/out.Bp2(end);
out.betat  = trapz(out.V,out.pres)/out.V(end)*2*4e-7*pi/out.bcentr^2;
out.beta   = out.betap*out.betat/(out.betap + out.betat);
out.betan  = abs(100*out.beta*out.amin*out.bcentr/(out.current/1e6));
out.ft     = (0.75*(1-out.h2.*(1-sqrt(1-out.h).*(1+0.5*out.h))./(out.h.^2)) + 0.25*(1-out.h2.*out.hh));
out.fc     = 1 - out.ft;
out.Li     = out.Bp2(end)*out.V(end)/((pi*4e-7)*out.current^2);
out.li     = 2*out.Li/(pi*4e-7*out.rcentr);
out.kappa95 = interp1(out.psin,out.kappa,0.95,'pchip');
out.deltau95 = interp1(out.psin,out.deltau,0.95,'pchip');
out.deltal95 = interp1(out.psin,out.deltal,0.95,'pchip');
out.delta95 = interp1(out.psin,out.delta,0.95,'pchip');
out.s95 = interp1(out.psin,out.s,0.95,'pchip');
out.q95 = interp1(out.psin,out.qpsi,0.95,'pchip');
out.jelite = -sign(out.current)*(out.rcentr*out.pprime.*(out.fpol/out.rcentr).^2.*out.B_2 + out.ffprim/Muo/out.rcentr);
out.jeliten = (max(out.jelite(out.psin>0.9)) + out.jelite(end))/2/(out.current/out.A(end));
%out.qstar = 2*pi*out.kappa95*out.amin^2*abs(out.bcentr)/out.rcentr/out.current/4e-7/pi;
%out.S_q4 = out.s./out.qpsi.^4;
%out.li3 = 2*out.Bp2(end)*out.V(end)/((pi*4e-7)^2*out.rcentr*out.current^2);

end