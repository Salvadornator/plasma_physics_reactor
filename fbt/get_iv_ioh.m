function params = get_iv_ioh(params)

%%
fprintf('       --- Calculating ohmic and vessel currents\n')
teq   = params.eqtime*1e-3;
dt    = params.dt;
Ip    = params.Ip;
Ia    = params.Ia;
Rp    = nanmin([params.vloop/Ip 1e-3]);
Lp    = params.mutuals.Lp;
Mpoh  = params.mutuals.Mpa(1);
Mpa   = params.mutuals.Mpa(2:17);
Mpv   = params.mutuals.Mpv;
Rvv   = diag(params.mutuals.Rv);
Mvv   = params.mutuals.Mvv;
Mvp   = params.mutuals.Mvp;
Mva   = params.mutuals.Mav.';
Mvoh  = Mva(:,1);
Mva   = params.mutuals.Mva(:,2:17);

Ip0   = params.Ip0;
Ioh0  = params.Ioh0;
Ia0   = params.Ia0;
Iv0   = params.Iv0;
Rp0   = nanmin([params.vloop0/Ip0 1e-3]);
Lp0   = params.Lp0;
Mpoh0 = params.Mpa0(1);
Mpa0  = params.Mpa0(2:17);
Mpv0  = params.Mpv0;
Mvp0  = params.Mvp0;

dIpdt   = (Ip - Ip0)/dt;
dIadt   = (Ia - Ia0)/dt;
dRpdt   = (Rp - Rp0)/dt;
dLpdt   = (Lp - Lp0)/dt;
dMvpdt  = (Mvp  - Mvp0)/dt;
dMpohdt = (Mpoh - Mpoh0)/dt;
dMpadt  = (Mpa  - Mpa0)/dt;
dMpvdt  = (Mpv  - Mpv0)/dt;

%%
if params.i_vessel
    L0 = -[Mpoh0 Mpv0; Mvoh Mvv];
    L1 = -[dMpohdt dMpvdt; zeros(160,1) zeros(160,160)];
    A  =  [dMpohdt dMpvdt; zeros(160,1) Rvv];
    B1 =  [Ip0*dRpdt + Rp0*dIpdt + 2*dIpdt*dLpdt + 2*dMpadt*dIadt; 2*dIpdt*dMvpdt];
    B2 =  [dIpdt*dRpdt; zeros(160,1)];
    S  =  [dMpadt*Ia0 + dLpdt*Ip0 + Lp0*dIpdt + Mpa0*dIadt + Ip0*Rp0; Ip0*dMvpdt + Mva*dIadt + Mvp0*dIpdt];
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,IohIv] = ode45(@(t,IohIv) get_vessel_curr(t,IohIv,L0,L1,A,B1,B2,S),[0 dt],[Ioh0; Iv0],opts);
    params.t_vv_oh =  t.' - t(end) + teq;
    params.Ioh = IohIv(:,1).';
    params.Iv  = IohIv(:,2:161).';
    t_aux = linspace(params.t_vv_oh(1),params.t_vv_oh(end),length(params.t_vv_oh));
    params.Ip_vv_oh = linspace(Ip0,Ip,length(params.t_vv_oh));
    params.Ip_vv_oh = interp1(t_aux,params.Ip_vv_oh,params.t_vv_oh,'pchip');
    params.Ia_vv_oh = zeros(length(Ia0),length(params.t_vv_oh));
    for ii = 1:length(Ia0)
        params.Ia_vv_oh(ii,:) = linspace(Ia0(ii),Ia(ii),length(params.t_vv_oh));
        params.Ia_vv_oh(ii,:) = interp1(t_aux,params.Ia_vv_oh(ii,:),params.t_vv_oh,'pchip');
    end
    dl = [params.dZv(1:44); params.dRv(45:80); params.dZv(81:124); params.dRv(125:160)];
    params.jv = params.Iv./repmat(dl,[1,length(params.t_vv_oh)]);
    params.Ive = params.mutuals.T_e_v*params.Iv;

    % Plotting
    c = 'kgrmb';
    figure(3)
    clf
    hold on
    for ii = 1:5
        plot(params.t_vv_oh*1e3,params.Ive(ii,:)/1e3,'--','color',c(ii))
    end
    plot(params.t_vv_oh*1e3,params.Ioh/1e3,'k','linewidth',3)
    plot(params.t_vv_oh*1e3,params.Ip_vv_oh/10e3,'r','linewidth',3)
    legend('I_{ve} ( #1 )','I_{ve} ( #2 )','I_{ve} ( #3 )','I_{ve} ( #4 )','I_{ve} ( #5 )','I_{OH}','I_P / 10','location','northwest')
    hold off
    xlabel('Time ( ms )')
    ylabel('Currents ( kA )')
    drawnow
elseif params.i_ioh
    params.t_vv_oh = linspace(teq-dt,teq,10001);
    dIohdt = -Mpoh\((Rp0 + dRpdt*(params.t_vv_oh - params.t_vv_oh(1))).*(Ip0 + dIpdt*(params.t_vv_oh - params.t_vv_oh(1))) + (Lp0 + dLpdt*(params.t_vv_oh - params.t_vv_oh(1)))*dIpdt + ...
                    (Ip0 + dIpdt*(params.t_vv_oh - params.t_vv_oh(1)))*dLpdt + (Mpa0*dIadt + dMpadt*dIadt*(params.t_vv_oh - params.t_vv_oh(1))) + (dMpadt*Ia0 + dMpadt*dIadt*(params.t_vv_oh - params.t_vv_oh(1))));
    params.Ioh = integration(params.t_vv_oh,dIohdt,params.Ioh0);
    
    % Plotting
    figure(3)
    clf
    plot(params.t_vv_oh*1e3,params.Ioh/1e3,'k','linewidth',2)
    xlabel('Time ( ms )')
    title('I_{OH} ( kA )')
    drawnow
end

end