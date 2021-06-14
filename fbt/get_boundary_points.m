function params = get_boundary_points(params)

%% Calculating boundary points
o = linspace(0,2*pi,params.Nb+1);
our = o; our(o>=pi/2) = 0;
oul = o; oul(o<pi/2) = 0; oul(o>=pi) = 0;
oll = o; oll(o<pi) = 0; oll(o>=3*pi/2) = 0;
olr = o; olr(o<3*pi/2) = 0;
Rup = params.Rmaxis + params.a*cos(o + params.deltau*sin(o) - params.sqrur*sin(2*our) - params.sqrul*sin(2*oul));
Rlo = params.Rmaxis + params.a*cos(o + params.deltal*sin(o) - params.sqrlr*sin(2*olr) - params.sqrll*sin(2*oll));
params.Rb = [Rup(o<=pi) Rlo(o>pi)];
params.Zb = params.Zmaxis + params.kappa*params.a.*sin(o);
params.Rb = params.Rb(1:end-1);
params.Zb = params.Zb(1:end-1);

end