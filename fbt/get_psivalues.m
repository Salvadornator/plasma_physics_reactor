function params = get_psivalues(params)

fprintf('       --- Calculating poloidal flux at prescribed points\n')

%% Computing coils poloidal flux at boundary points
params.fbt.psibp = interp2(params.rmesh,params.zmesh,params.psip,[params.Rb params.Rxp params.Rsf params.Rbfix],[params.Zb params.Zxp params.Zsf params.Zbfix],'cubic').';
params.fbt.Gbc = zeros(length([params.Rb params.Rxp params.Rsf params.Rbfix]),size(params.coils,1));
for icoil = 1:params.Nc_E
    eval(['params.fbt.Gbc(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,[params.Rb params.Rxp params.Rsf params.Rbfix],[params.Zb params.Zxp params.Zsf params.Zbfix],''cubic'').'';'])
end
for icoil = 1:params.Nc_F
    eval(['params.fbt.Gbc(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,[params.Rb params.Rxp params.Rsf params.Rbfix],[params.Zb params.Zxp params.Zsf params.Zbfix],''cubic'').'';'])
end
for icoil = 1:params.Nc_D
    eval(['params.fbt.Gbc(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,[params.Rb params.Rxp params.Rsf params.Rbfix],[params.Zb params.Zxp params.Zsf params.Zbfix],''cubic'').'';'])
end
params.fbt.Gboh = interp2(params.greens.r,params.greens.z,-params.greens.OH.G/2/pi,[params.Rb params.Rxp params.Rsf params.Rbfix],[params.Zb params.Zxp params.Zsf params.Zbfix],'cubic').';
RRB = [params.Rb params.Rxp params.Rsf params.Rbfix];
ZZB = [params.Zb params.Zxp params.Zsf params.Zbfix];
params.fbt.Gbv = zeros(size(params.vv,1),length(RRB));
for ibp = 1:length(ZZB)
    params.fbt.Gbv(:,ibp) = green_em(RRB(ibp),ZZB(ibp),params.vv(:,1),params.vv(:,2),'mutual');
end
params.fbt.Gbv = -params.fbt.Gbv.'/2/pi;

%% Computing coils poloidal flux at fixed boundary points
if ~isempty(params.Rbfix)
    params.fbt.psibfixp = interp2(params.rmesh,params.zmesh,params.psip,params.Rbfix,params.Zbfix,'cubic').';
    params.fbt.Gbfixc = zeros(length(params.Rbfix),size(params.coils,1));
    for icoil = 1:params.Nc_E
        eval(['params.fbt.Gbfixc(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,params.Rbfix,params.Zbfix,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_F
        eval(['params.fbt.Gbfixc(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,params.Rbfix,params.Zbfix,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_D
        eval(['params.fbt.Gbfixc(:,' num2str(icoil+params.Nc_E+params.Nc_D) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,params.Rbfix,params.Zbfix,''cubic'').'';'])
    end
    params.fbt.Gbfixoh = interp2(params.greens.r,params.greens.z,-params.greens.OH.G/2/pi, params.Rbfix,params.Zbfix,'cubic').';
    params.fbt.Gbfixv = zeros(size(params.vv,1),length(params.Rbfix));
    for ibp = 1:size(params.vv,1)
        params.fbt.Gbfixv(ibp,:) = green_em(params.vv(ibp,1),params.vv(ibp,2),params.Rbfix,params.Zbfix,'mutual');
    end
    params.fbt.Gbfixv = -params.fbt.Gbfixv.'/2/pi;
else
    params.fbt.psibfixp = [];
    params.fbt.Gbfixc = [];
    params.fbt.Gbfixoh = 0;
    params.fbt.Gbfixv = zeros(1,160);
end

%% Computing coils poloidal flux at magnetic axis
if isfield(params,'geqdsk')
    rmaxis = params.geqdsk.rmaxis;
    zmaxis = params.geqdsk.zmaxis;
else
    rmaxis = params.Rmaxis;
    zmaxis = params.Zmaxis;
end
params.fbt.Gcmaxis = zeros(1,size(params.coils,1));
for icoil = 1:params.Nc_E
    eval(['params.fbt.Gcmaxis(1,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,rmaxis,zmaxis,''cubic'').'';'])
end
for icoil = 1:params.Nc_F
    eval(['params.fbt.Gcmaxis(1,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,rmaxis,zmaxis,''cubic'').'';'])
end
for icoil = 1:params.Nc_D
    eval(['params.fbt.Gcmaxis(1,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,rmaxis,zmaxis,''cubic'').'';'])
end

%% Computing coils poloidal flux at x-points
if ~isempty(params.Rxp)
    params.fbt.psixpp = interp2(params.greens.r,params.greens.z,params.psip,params.Rxp,params.Zxp,'cubic').';
    params.fbt.Gxpc = zeros(length(params.Rxp),size(params.coils,1));
    for icoil = 1:params.Nc_E
        eval(['params.fbt.Gxpc(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_F
        eval(['params.fbt.Gxpc(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_D
        eval(['params.fbt.Gxpc(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    if strcmp(params.plasmaid,'snowflake')
        params.fbt.Gxpc = params.fbt.Gxpc(1,:);
        params.fbt.psixpp = params.fbt.psixpp(1);
    end
    params.fbt.Gxpoh = interp2(params.greens.r,params.greens.z,-params.greens.OH.G/2/pi, params.Rxp,params.Zxp,'cubic').';
    params.fbt.Gxpv = zeros(size(params.vv,1),length(params.Rxp));
    for ibp = 1:length(params.Rxp)
        params.fbt.Gxpv(:,ibp) = green_em(params.Rxp(ibp),params.Zxp(ibp),params.vv(:,1),params.vv(:,2),'mutual');
    end
    params.fbt.Gxpv = -params.fbt.Gxpv.'/2/pi;
else
    params.fbt.psixpp = [];
    params.fbt.Gxpc = [];
    params.fbt.Gxpoh = 0;
    params.fbt.Gxpv = zeros(1,160);
end

%% Computing coils poloidal flux at second-order null
if ~isempty(params.Rsf)
    params.fbt.psisfp = interp2(params.greens.r,params.greens.z,params.psip,params.Rsf,params.Zsf,'cubic').';
    params.fbt.Gsfc = zeros(length(params.Rsf),size(params.coils,1));
    for icoil = 1:params.Nc_E
        eval(['params.fbt.Gsfc(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,params.Rsf,params.Zsf,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_F
        eval(['params.fbt.Gsfc(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,params.Rsf,params.Zsf,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_D
        eval(['params.fbt.Gsfc(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,params.Rsf,params.Zsf,''cubic'').'';'])
    end
    params.fbt.Gsfoh = interp2(params.greens.r,params.greens.z,-params.greens.OH.G/2/pi, params.Rsf,params.Zsf,'cubic').';
    params.fbt.Gsfv = zeros(size(params.vv,1),length(params.Rsf));
    for ibp = 1:length(params.Rsf)
        params.fbt.Gsfv(:,ibp) = green_em(params.Rsf(ibp),params.Zsf(ibp),params.vv(:,1),params.vv(:,2),'mutual');
    end
    params.fbt.Gsfv = -params.fbt.Gsfv.'/2/pi;
else
    params.fbt.psisfp = [];
    params.fbt.Gsfc = [];
    params.fbt.Gsfoh = 0;
    params.fbt.Gsfv = zeros(1,160);
end

%% Computing coils poloidal flux at strike points
if ~isempty(params.Rsp)
    params.fbt.psispp = interp2(params.greens.r,params.greens.z,params.psip,params.Rsp,params.Zsp,'cubic').';
    params.fbt.Gspc = zeros(length(params.Rsp),size(params.coils,1));
    for icoil = 1:params.Nc_E
        eval(['params.fbt.Gspc(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,-params.greens.E' num2str(icoil) '.G/2/pi,params.Rsp,params.Zsp,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_F
        eval(['params.fbt.Gspc(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,-params.greens.F' num2str(icoil) '.G/2/pi,params.Rsp,params.Zsp,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_D
        eval(['params.fbt.Gspc(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,-params.greens.D' num2str(icoil) '.G/2/pi,params.Rsp,params.Zsp,''cubic'').'';'])
    end
    params.fbt.Gspoh = interp2(params.greens.r,params.greens.z,-params.greens.OH.G/2/pi, params.Rsp,params.Zsp,'cubic').';
    params.fbt.Gspv = zeros(size(params.vv,1),length(params.Rsp));
    for ibp = 1:length(params.Rsp)
        params.fbt.Gspv(:,ibp) = green_em(params.Rsp(ibp),params.Zsp(ibp),params.vv(:,1),params.vv(:,2),'mutual');
    end
    params.fbt.Gspv = -params.fbt.Gspv.'/2/pi;
else
    params.fbt.psispp = [];
    params.fbt.Gspc = [];
    params.fbt.Gspoh = 0;
    params.fbt.Gspv = zeros(1,160);
end

%% Computing BR and BZ at x-points due to the plasma
if ~isempty(params.Rxp)
    params.fbt.brxpp = interp2(params.greens.r,params.greens.z,params.brp,params.Rxp,params.Zxp,'cubic').';
    params.fbt.bzxpp = interp2(params.greens.r,params.greens.z,params.bzp,params.Rxp,params.Zxp,'cubic').';
    params.fbt.brcxp = zeros(length(params.Rxp),size(params.coils,1));
    params.fbt.bzcxp = zeros(length(params.Rxp),size(params.coils,1));
    for icoil = 1:10
        eval(['params.fbt.brcxp(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.BR,params.Rxp,params.Zxp,''cubic'').'';'])
        eval(['params.fbt.bzcxp(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.BZ,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    for icoil = 1:4
        eval(['params.fbt.brcxp(:,' num2str(icoil+10) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.BR,params.Rxp,params.Zxp,''cubic'').'';'])
        eval(['params.fbt.bzcxp(:,' num2str(icoil+10) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.BZ,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    for icoil = 1:2
        eval(['params.fbt.brcxp(:,' num2str(icoil+14) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.BR,params.Rxp,params.Zxp,''cubic'').'';'])
        eval(['params.fbt.bzcxp(:,' num2str(icoil+14) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.BZ,params.Rxp,params.Zxp,''cubic'').'';'])
    end
    params.fbt.brxpoh = interp2(params.greens.r,params.greens.z,params.greens.OH.BR, params.Rxp,params.Zxp,'cubic').';
    params.fbt.bzxpoh = interp2(params.greens.r,params.greens.z,params.greens.OH.BZ, params.Rxp,params.Zxp,'cubic').';
    params.fbt.brxpv = zeros(size(params.vv,1),length(params.Rxp));
    params.fbt.bzxpv = zeros(size(params.vv,1),length(params.Rxp));
    for ibp = 1:size(params.vv,1)
        params.fbt.brxpv(ibp,:) = green_em(params.vv(ibp,1),params.vv(ibp,2),params.Rxp,params.Zxp,'br');
        params.fbt.bzxpv(ibp,:) = green_em(params.vv(ibp,1),params.vv(ibp,2),params.Rxp,params.Zxp,'bz');
    end
    params.fbt.brxpv = params.fbt.brxpv.';
    params.fbt.bzxpv = params.fbt.bzxpv.';
else
    params.fbt.brxpp = [];
    params.fbt.bzxpp = [];
    params.fbt.brcxp = [];
    params.fbt.bzcxp = [];
    params.fbt.brxpoh = 0;
    params.fbt.bzxpoh = 0;
    params.fbt.brxpv = zeros(1,160);
    params.fbt.bzxpv = zeros(1,160);
end

%% Computing BR, BZ, dBRdR, dBRdz, dBZdR and dBZdZ at second-order null due to the plasma
if ~isempty(params.Rsf)
    params.fbt.brsfp = interp2(params.greens.r,params.greens.z,params.greens.brp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.bzsfp = interp2(params.greens.r,params.greens.z,params.greens.bzp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.brcsf = zeros(length(params.Rsf),size(params.coils,1));
    params.fbt.bzcsf = zeros(length(params.Rsf),size(params.coils,1));
    params.fbt.dbrdrsfp = interp2(params.greens.r,params.greens.z,params.greens.dbrdrp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.dbrdzsfp = interp2(params.greens.r,params.greens.z,params.greens.dbrdzp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.dbzdrsfp = interp2(params.greens.r,params.greens.z,params.greens.dbzdrp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.dbzdzsfp = interp2(params.greens.r,params.greens.z,params.greens.dbzdzp,params.Rsf,params.Zsf,'cubic').';
    params.fbt.dbrdrcsf = zeros(length(params.Rsf),size(params.coils,1));
    params.fbt.dbrdzcsf = zeros(length(params.Rsf),size(params.coils,1));
    params.fbt.dbzdrcsf = zeros(length(params.Rsf),size(params.coils,1));
    params.fbt.dbzdzcsf = zeros(length(params.Rsf),size(params.coils,1));
    for icoil = 1:params.Nc_E
        eval(['params.fbt.brcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.BR.,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.bzcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.BZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdrcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.dBRdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdrcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.dBZdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdzcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.dBRdZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdzcsf(:,' num2str(icoil) ') = interp2(params.greens.r,params.greens.z,params.greens.E' num2str(icoil) '.dBZdZ,params.Rsf,params.Zsf,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_F
        eval(['params.fbt.brcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.BR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.bzcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.BZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdrcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.dBRdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdrcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.dBZdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdzcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.dBRdZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdzcsf(:,' num2str(icoil+params.Nc_E) ') = interp2(params.greens.r,params.greens.z,params.greens.F' num2str(icoil) '.dBZdZ,params.Rsf,params.Zsf,''cubic'').'';'])
    end
    for icoil = 1:params.Nc_D
        eval(['params.fbt.brcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.BR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.bzcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.BZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdrcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.dBRdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdrcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.dBZdR,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbrdzcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.dBRdZ,params.Rsf,params.Zsf,''cubic'').'';'])
        eval(['params.fbt.dbzdzcsf(:,' num2str(icoil+params.Nc_E+params.Nc_F) ') = interp2(params.greens.r,params.greens.z,params.greens.D' num2str(icoil) '.dBZdZ,params.Rsf,params.Zsf,''cubic'').'';'])
    end
else
    params.fbt.brsfp = [];
    params.fbt.bzsfp = [];
    params.fbt.brcsf = [];
    params.fbt.bzcsf = [];
    params.fbt.dbrdrsfp = [];
    params.fbt.dbzdrsfp = [];
    params.fbt.dbrdzsfp = [];
    params.fbt.dbzdzsfp = [];
    params.fbt.dbrdrcsf = [];
    params.fbt.dbzdrcsf = [];
    params.fbt.dbrdzcsf = [];
    params.fbt.dbzdzcsf = [];
end

end