function params = j0_guess(params,iplot)

%% Defining grid
fprintf('   --- Plasma current density: initial guess\n')
params.jcur.rg = linspace(params.rmesh(1),params.rmesh(end),33).';
params.jcur.zg = linspace(params.zmesh(1),params.zmesh(end),33).'; % params.gridsize
[Rg,~] = meshgrid(params.jcur.rg,params.jcur.zg);

%% Calculating boundary points
xx = 1;
o = linspace(0,2*pi,51);
our = o; our(o>=pi/2) = 0;
oul = o; oul(o<pi/2) = 0; oul(o>=pi) = 0;
oll = o; oll(o<pi) = 0; oll(o>=3*pi/2) = 0;
olr = o; olr(o<3*pi/2) = 0;
Rup = params.Rmaxis + params.a*xx.*cos(o + (params.deltau*xx.^2)*sin(o) - (params.sqrur*xx.^2)*sin(2*our) - (params.sqrul*xx.^2)*sin(2*oul));
Rlo = params.Rmaxis + params.a*xx.*cos(o + (params.deltal*xx.^2)*sin(o) - (params.sqrlr*xx.^2)*sin(2*olr) - (params.sqrll*xx.^2)*sin(2*oll));
params.Rbb = [Rup(o<=pi) Rlo(o>pi)]; params.Rbb = params.Rbb(1:end-1);
params.Zbb = params.Zmaxis + (1 + (params.kappa-1)*xx.^2)*params.a.*xx.*sin(o);  params.Zbb = params.Zbb(1:end-1);

rup = params.Rmaxis + params.a*0.75.*cos(o + (params.deltau*0.75.^2)*sin(o) - (params.sqrur*0.75.^2)*sin(2*our) - (params.sqrul*0.75.^2)*sin(2*oul));
rlo = params.Rmaxis + params.a*0.75.*cos(o + (params.deltal*0.75.^2)*sin(o) - (params.sqrlr*0.75.^2)*sin(2*olr) - (params.sqrll*0.75.^2)*sin(2*oll));
params.jcur.rbbj = [rup(o<=pi) rlo(o>pi)];
params.jcur.zbbj = params.Zmaxis + (1 + (params.kappa-1)*0.75.^2)*params.a.*0.75.*sin(o);

%% Computing J distribution
J = zeros(size(Rg));
for ir = 1:length(params.jcur.rg)
    for iz = 1:length(params.jcur.zg)
        if inpolygon(params.jcur.rg(ir),params.jcur.zg(iz),params.jcur.rbbj,params.jcur.zbbj)
            J(iz,ir)  = 1;
        end
    end
end

%% Normalizing J to the plasma current
dR = mean(diff(params.jcur.rg));
dZ = mean(diff(params.jcur.zg));
params.jcur.jplasma = J/(sum(J(:)*dR*dZ))*params.Ip;
%[R,Z] = meshgrid(params.rmesh,params.zmesh);
%ind = inpolygon(R,Z,params.Rbb,params.Zbb);
%params.jplasma = interp2(params.jcur.rg,params.jcur.zg,params.jcur.jplasma,R,Z);
%params.jplasma(~ind) = 0;

%% Plotting J distribution
if iplot
    figure(1)
    grid on
    clf
    contourf(params.jcur.rg,params.jcur.zg,params.jcur.jplasma/1e6,21)
    hold on
    plot(params.Rbb,params.Zbb,'r','linewidth',3)
    plot(params.Rmaxis,params.Zmaxis,'w+','linewidth',3,'markersize',10)
    plot(params.Rxp,params.Zxp,'wx','linewidth',3,'markersize',10)
    plot_tcabr('fig',1)
    hold off
    title('Plasma Current Density Distribution ( MA / m^2 )')
    caxis([0 ceil(max(params.jcur.jplasma(:)/1e6))-1])
    colorbar
    axis([0.35 0.9 -0.3001 0.3001])
    drawnow
end

end