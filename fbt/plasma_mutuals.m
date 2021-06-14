function params = plasma_mutuals(params)

fprintf('       --- Calculating mutual inductances between plasma and machine\n')

% J-distribution
dr = mean(diff(params.jcur.rg));
dz = mean(diff(params.jcur.zg));
[R,Z] = meshgrid(params.jcur.rg,params.jcur.zg);
Ix = params.jcur.jplasma(:)*dr*dz;

% Plasma self-inductance
Rx = R(:);
Zx = Z(:);
Mpp = zeros(length(Rx),length(Rx));
for ix = 1:length(Rx)
    Mpp(ix,:) = green_em(Rx(ix),Zx(ix),Rx,Zx,'mutual');
end
Lp = (Ix.'*Mpp*Ix)/params.Ip^2;

% Mutual between coils and plasma
coilnames = {'OH','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','F1','F2','F3','F4','D1','D2','D3','D4'};
Map = zeros(length(coilnames),length(Rx));
for ic = 1:length(coilnames)
    eval(['Map(ic,:) = interp2(params.greens.r,params.greens.z,params.greens.' coilnames{ic} '.G,Rx,Zx,''cubic'');'])
end
Mpa = Map.';
Mpa = (Ix.'*Mpa)/params.Ip;
Map = Mpa.';

% Mutual between vessel and plasma
Mvp = zeros(size(params.vv,1),length(Rx));
for ivessel = 1:size(params.vv,1)
    Mvp(ivessel,:) = green_em(params.vv(ivessel,1),params.vv(ivessel,2),Rx,Zx,'mutual').';
end
Mpv = Mvp.';
Mpv = (Ix.'*Mpv)/params.Ip;
Mvp = Mpv.';

%% Outputs
params.mutuals.Map = Map;
params.mutuals.Mpa = Mpa;
params.mutuals.Mpv = Mpv;
params.mutuals.Mvp = Mvp;
params.mutuals.Lp  = Lp;
params.mutuals.Rp  = nanmin([params.vloop/params.Ip 1e-3]);

end