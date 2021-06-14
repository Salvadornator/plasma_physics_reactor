function params = new_jplasma(params,iplot)

fprintf('       --- Updating plasma current density\n')

  %% Updating the J-profile
  if strcmp(params.imode,'hmode')
      struct = pedestal_profiles(params);
      params.u_hat = struct.u_hat;
      params.omega_hat = struct.omega_hat;
      params.omega_exb = struct.omega_exb;
      params.utor = struct.utor;
      params.upol = struct.upol;
      params.vloop = struct.Vl;
      params.Pohm = struct.Pohm;
      params.Prad_H = struct.Prad_H;
      params.Prad_C = struct.Prad_C;
      params.Prad_Fe = struct.Prad_Fe;
      params.Prad_brem_H = struct.Prad_brem_H;
      params.Prad_brem_C = struct.Prad_brem_C;
      params.Prad_brem_Fe = struct.Prad_brem_Fe;
      params.Prad = struct.Prad;
      p = struct.ptot;
      pprime  = struct.pprime;
      ffprime = struct.ffprim;
      f = sign(params.Bcentr)*sqrt((params.Rcentr*params.Bcentr)^2 + cumtrapz(struct.psin,2*ffprime*(struct.sibry-struct.simag)));
      [R,Z]    = meshgrid(params.geqdsk.rg,params.geqdsk.zg);
      psirzn   = params.geqdsk.psirzn;
      ind      = inpolygon(R,Z,params.geqdsk.rbbbs,params.geqdsk.zbbbs);
      psirzn(~ind) = 0;
      pprim   = reshape(interp1(params.geqdsk.psin,pprime,psirzn(:)),params.gridsize,params.gridsize);
      ffprim  = reshape(interp1(params.geqdsk.psin,ffprime,psirzn(:)),params.gridsize,params.gridsize);
      jplasma = -sign(struct.current)*(R.*pprim + ffprim./R/(pi*4e-7)); jplasma(~ind) = 0;
      %icur = sum(jplasma(:)*mean(diff(params.geqdsk.rg))*mean(diff(params.geqdsk.zg)))
      jtor    = -sign(struct.current)*(struct.rcentr*pprime + ffprime/struct.rcentr/(pi*4e-7).*struct.R02_R2);
      jpol    = -sign(struct.current)*ffprime/struct.rcentr/(pi*4e-7).*struct.Bp2./(f/struct.rcentr).^2;
      q       = -sign(params.Ip)*f.*params.geqdsk.Vprime.*params.geqdsk.R2_1/4/pi^2;

      % Calculating line average density
      zline = linspace(min(params.geqdsk.zg),max(params.geqdsk.zg),10001)';
      rline = params.Rcentr*ones(size(zline));
      neavg   = reshape(interp1(struct.psin,struct.ne,psirzn(:)),params.gridsize,params.gridsize);
      neavg(~ind) = 0;
      neavg   = interp2(params.geqdsk.rg,params.geqdsk.zg,neavg,rline,zline);
      Lne = inpolygon(rline,zline,params.geqdsk.rbbbs,params.geqdsk.zbbbs);
      params.Lne = max(zline(Lne)) - min(zline(Lne));
      params.neavg = trapz(zline,neavg)/params.Lne;
      params.nG    = params.geqdsk.nG;
      params.P_LH = 0.072e6*((params.neavg/1e20)^0.7).*((abs(params.geqdsk.bcentr))^0.7).*((params.geqdsk.Sperp(end))^0.9).*(mean(struct.Zeff/2)^0.7).*(0.1*(params.geqdsk.rcentr/params.geqdsk.amin)./(1-sqrt(2/(1+params.geqdsk.rcentr/params.geqdsk.amin)))).^0.5;
      %params.P_LH = 0.042e6*((params.neavg/1e20)^0.73).*((abs(params.geqdsk.bcentr))^0.74).*((params.geqdsk.Sperp(end))^0.98);
      params.fce = params.geqdsk.fce;
      params.rho_ce = sqrt(1.6e-19*struct.Te/9.11e-31)./params.fce/2/pi;
      params.fci = params.geqdsk.fci;
      params.rho_ci = sqrt(1.6e-19*struct.Ti/1.67e-27)./params.fci/2/pi;
      params.rho_star = interp1(params.geqdsk.rho_tor,params.rho_ci/params.geqdsk.amin,0.5,'pchip');
      
      % Alpha particles
      sv = sigmav(struct.Ti);
      p_fus = 0.25*1.602e-19*17.6e6*struct.ni.^2.*sv;
      params.Pfus = trapz(struct.V,p_fus);
      
      fprintf('       ---   Pohm: %3d kW\n',round(struct.Pohm/1e3))
      fprintf('       ---   Paux: %3d kW\n',round(params.Paux/1e3))
      fprintf('       ---   Prad: %3d kW\n',round(struct.Prad/1e3))
      fprintf('       ---   Psol: %3d kW\n',round((struct.Pohm + params.Paux + params.Pfus/5 - struct.Prad)/1e3))
      fprintf('       ---   P_LH: %3d kW\n',round(params.P_LH/1e3))
      fprintf('       ---   Pfus: %3d kW\n',round(params.Pfus/1e3))
      fprintf('       ---   Qfus: %4.2f\n',params.Pfus/(struct.Pohm + params.Paux))
      fprintf('       ---  betaN: %4.2f\n',params.geqdsk.betan)
      fprintf('       ---  neavg: %4.2fx10^20 m^-3\n',params.neavg/1e20)
      fprintf('       ---     nG: %4.2fx10^20 m^-3\n',params.nG/1e20)
      fprintf('       ---  rhost: %4.2f %%\n',params.rho_star*100)
      fprintf('       ---  Vloop: %4.2f Volts\n',params.vloop)
  else
      [p,pprime,f,ffprime,jplasma] = gs_global_params2(params);
  end
  params.geqdsk.pprime  = pprime;
  params.geqdsk.ffprim = ffprime;
  params.geqdsk.fpol = f;
  params.geqdsk.pres = p;
  params.geqdsk.qpsi = q;
  params.geqdsk.jpol = jpol;
  params.geqdsk.jtor = jtor;
  params.geqdsk.jpar = struct.jpar;
  params.jplasma = jplasma;
  [Rj,Zj] = meshgrid(params.jcur.rg,params.jcur.zg);
  params.jcur.jplasma = interp2(params.rmesh,params.zmesh,params.jplasma,Rj,Zj);
  in = inpolygon(params.jcur.rg,params.jcur.zg,params.geqdsk.rbbbs,params.geqdsk.zbbbs);
  params.jcur.jplasma(~in) = 0;
  writeg(params.geqdsk,'geqdsk_out')
  params.geqdsk = readg('./geqdsk_out');
  params.tau_E = get_tau_e(params.geqdsk.current/1e6,params.geqdsk.bcentr,(params.Pohm + params.Paux + params.Pfus/5)/1e6,params.neavg/1e19,1,params.geqdsk.rcentr,params.geqdsk.amin,params.geqdsk.kappa95);
  fprintf('       ---  tau_E: %4.2f ms\n',params.tau_E*1e3)
  fprintf('       --- Wp (GEQDSK file): %4.2f kJ\n',params.geqdsk.wp/1e3)
  fprintf('       --- Wp (scaling law): %4.2f kJ\n',params.tau_E*(params.Pohm + params.Paux + params.Pfus/5 - params.Prad)/1e3)
  
  % Plotting
      figure(4)
      clf
      hs(1) = subplot(3,2,1);
      plot(params.geqdsk.psin,params.geqdsk.pres/1e3,'k','linewidth',2)
      xlabel('\psi_N')
      legend('P ( kPa )','location','southwest')
      hs(2) = subplot(3,2,3);
      plot(params.geqdsk.psin,params.geqdsk.fpol,'k','linewidth',2)
      legend('F ( T\cdotm )','location','northwest')
      xlabel('\psi_N')
      hs(3) = subplot(3,2,2);
      plot(params.geqdsk.psin,params.geqdsk.rcentr*params.geqdsk.pprime/1e6,'k','linewidth',2)
      xlabel('\psi_N')
      legend('R_0 P'' ( MA / m^2 )','location','southwest')
      hs(4) = subplot(3,2,4);
      plot(params.geqdsk.psin,params.geqdsk.ffprim/params.geqdsk.rcentr/(pi*4e-7)/1e6,'k','linewidth',2)
      legend('FF'' / \mu_0 R_0 ( MA / m^2 )','location','southeast')
      xlabel('\psi_N')
      hs(5) = subplot(3,2,5);
      plot(params.geqdsk.psin,params.geqdsk.qpsi,'k','linewidth',2)
      hold on
      plot(params.geqdsk.psin,params.geqdsk.s,'color',[0 0.5 0],'linewidth',2)
      hold off
      legend('q','s','location','northwest')
      xlabel('\psi_N')
      axis([0 1 0 1.5*params.geqdsk.q95])
      hs(6) = subplot(3,2,6);
      plot(params.geqdsk.psin,params.geqdsk.jtor/1e6,'k','linewidth',2)
      legend('J_{tor} ( MA / m^2 )')
      xlabel('\psi_N')
      linkaxes(hs,'x')

  %% Plotting J distribution
if iplot
    figure(1)
    grid on
    clf
    [~,h] = contourf(params.rmesh,params.zmesh,params.jplasma,51);
    set(h,'LineColor','none')
    hold on
    plot(params.Rbb,params.Zbb,'r','linewidth',3)
    plot(params.fbt.rE,params.fbt.zE,'w+','linewidth',3,'markersize',10)
    plot(params.Rxp,params.Zxp,'wx','linewidth',3,'markersize',10)
    plot_tcabr
    hold off
    title('Plasma Current Density Distribution ( MA / m^2 )')
    %caxis([0 ceil(max(params.jplasma(:)/1e6))-1])
    drawnow
end

end