function p = getprofs(p)

pathname = p.kinprofs;
p.ptot = load(fullfile(pathname,'profile_p'));    p.ptot = p.ptot(:,2)*1e3;
p.ne   = load(fullfile(pathname,'profile_ne'));   p.ne = p.ne(:,2)*1e20;
p.ni   = load(fullfile(pathname,'profile_ni'));   p.ni = p.ni(:,2)*1e20;
p.te   = load(fullfile(pathname,'profile_te'));   p.te = p.te(:,2)*1e3;
p.ti   = load(fullfile(pathname,'profile_ti'));   p.ti = p.ti(:,2)*1e3;
p.nc   = load(fullfile(pathname,'profile_nc'));   p.nc = p.nc(:,2)*1e20;
p.nFe  = load(fullfile(pathname,'profile_nFe'));  p.nFe = p.nFe(:,2)*1e20;
p.zeff = load(fullfile(pathname,'profile_zeff')); p.zeff = p.zeff(:,2);
p.fz   = load(fullfile(pathname,'profile_fz'));   p.fz = p.fz(:,2);

end