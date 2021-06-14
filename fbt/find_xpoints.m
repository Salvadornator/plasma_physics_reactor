function output = find_xpoints(params)

fprintf('       --- Finding magnetic axis and x-points\n')
r = params.rmesh;
z = params.zmesh;
psipol = params.psi;
Ip = params.Ip;
Ro = params.Rmaxis;
Zo = params.Zmaxis;
a = params.a;

%% Loading generic GEQDSK
if isfield(params,'geqdsk')
    g = params.geqdsk;
else
    if params.Rcentr == 0.62
        g = readg('/home/canal/matlab/equil_reconst/fbt/geqdsk_in');
    else
        g = readg('/home/fmsalvador/matlab/plasma_physics_project/reactor/scenarios/geqdsk_reactor_in');
    end
end
if isempty(g)
    output = [];
    return
end
psin = linspace(0,1,length(r));
g.fpol = interp1(g.psin,g.fpol,psin).';
g.Ipol = interp1(g.psin,g.Ipol,psin).';
g.pres = interp1(g.psin,g.pres,psin).';
g.ffprim = interp1(g.psin,g.ffprim,psin).';
g.pprime = interp1(g.psin,g.pprime,psin).';
g.qpsi = interp1(g.psin,g.qpsi,psin).';
g.psirz = psipol;
g.bcentr = params.Bcentr;
g.rcentr = params.Rcentr;
g.rg = r.';
g.zg = z.';
g.nw = length(r);
g.nh = length(z);
g.rleft = r(1);
g.rdim = r(end) - r(1);
g.zmid = 0;
g.zdim = z(end) - z(1);
g.bcentr = params.Bcentr;
g.rcentr = params.Rcentr;
g.current = Ip;
g.rhovn = sqrt(psin);
g.psin = psin;
g.nw = length(psin);
g.nh = length(psin);
g.rlim = interp1(1:length(params.limiter(:,1)),params.limiter(:,1),linspace(1,length(params.limiter(:,1)),86));
g.zlim = interp1(1:length(params.limiter(:,2)),params.limiter(:,2),linspace(1,length(params.limiter(:,2)),86));
g.limitr = length(g.rlim);
if isfield(params.fbt,'fE') && ~isempty(params.Rxp) && params.iter > 3
    g.simag = params.fbt.fE;
else
    [R,Z] = meshgrid(params.rmesh,params.zmesh);
    ind = inpolygon(R,Z,params.Rbb,params.Zbb);
    psi_in = g.psirz;
    psi_in(~ind) = nan;
    [g.simag,indmax] = nanmin(psi_in(:));
    [indz,indr]=ind2sub(size(g.psirz),indmax);
    g.rmaxis = g.rg(indr);
    g.zmaxis = g.zg(indz);
end

if isfield(params.fbt,'fS') && ~isempty(params.Rxp) && params.iter > 3
    g.sibry = min(params.fbt.fS);
else
    g.sibry = min(interp2(g.rg,g.zg,g.psirz,interp1(1:length(params.limiter(:,1)),params.limiter(:,1),linspace(1,length(params.limiter(:,1)),1001)),interp1(1:length(params.limiter(:,2)),params.limiter(:,2),linspace(1,length(params.limiter(:,2)),1001))));
end
g.psirzn = (g.psirz - g.simag)/(g.sibry - g.simag);

%% Preparing the GEQDSK file into a PsiTbx function
writeg(g,'geqdsk_out')
g = readg('geqdsk_out');

if isempty(g)
    fprintf('       --- Error reading GEQDSK file: line 61 in find_xpoints.m\n\n')
    output = [];
    return
end
out = read_eqdsk('geqdsk_out',5);
delete('geqdsk_out')
psitbx = eqdsk2psitbx(out);

%% Finding the magnetic axis and x-points
if strcmp(params.gridtype,'large') && params.gridsize <= 129
    disp('       --- Increasing resolution of the grid: x2')
    [re,ze,~,rs,zs,~] = findcpts(psitbx,'resolution','x2');
else
    [re,ze,~,rs,zs,~] = findcpts(psitbx);
end

%% Selecting the points inside the vacuum vessel
in = inpolygon(re,ze,params.limiter(:,1),params.limiter(:,2));
rE = re(in); zE = ze(in); fE = interp2(g.rg,g.zg,g.psirz,rE,zE,'cubic');
in = inpolygon(rs,zs,params.limiter(:,1),params.limiter(:,2));
rS = rs(in); zS = zs(in); fS = interp2(g.rg,g.zg,g.psirz,rS,zS,'cubic');

%% Checking for expurious x-points/magnetic axes inside/outside the plasma
oo = linspace(0,2*pi,101);
Rc = Ro + a/2*cos(oo);
Zc = Zo + a/2*sin(oo);
in = inpolygon(rE,zE,Rc,Zc);
rE = rE(in); zE = zE(in); fE = fE(in);
in = inpolygon(rS,zS,Rc,Zc);
rS = rS(~in); zS = zS(~in); fS = fS(~in);

%% Minimum psi along the vacuum vessel wall
rwall = interp1(1:length(params.limiter(:,1)),params.limiter(:,1),1:.0001:length(params.limiter(:,1)));
zwall = interp1(1:length(params.limiter(:,1)),params.limiter(:,2),1:.0001:length(params.limiter(:,1)));
psi_wall = interp2(g.rg,g.zg,g.psirz,rwall,zwall,'cubic');
psi_wall = min(psi_wall);

%% Detecting if plasma is limited or diverted
psiminimum = min(g.psirz(:));
if isempty(rS) % No x-points inside the vacuum vessel
    psi_sep = psi_wall;
    Cw = contourc(g.rg,g.zg,g.psirz-psiminimum,([1 1]*psi_sep-psiminimum)*0.999);
    fcc = find_closed_contour(Cw,params);
    rbbbs = fcc.rb;
    zbbbs = fcc.zb;
    plasma = 'limited';
else % If there are x-points in the vacuum vessel, is the plasma diverted or limited?
    psi_sep = min(fS); % If plasma is diverted, the x-point with max flux determines the separarix
    CQ = contourc(g.rg,g.zg,g.psirz-psiminimum,([1 1]*psi_sep-psiminimum)*0.999);
    fcc = find_closed_contour(CQ,params);
    rbbbs = fcc.rb;
    zbbbs = fcc.zb;
    plasma = 'diverted';
    if isempty(fcc.rb) % If no closed flux contour is found for the x-points, then x-points are not active and the plasma is limited
        plasma = 'limited';
        if params.Rcentr == 0.62
            %params_fake.limiter = [[0.395;0.395;0.85;0.933;0.933;0.85;0.395] [-0.372;0.372;0.372;0.134;-0.134;-0.372;-0.372]];
            params_fake.limiter = [[params.limiter(:,1);params.limiter(1,1)] [params.limiter(:,2);params.limiter(1,2)]];
        else
            params_fake.limiter = [[params.limiter(:,1);params.limiter(1,1)] [params.limiter(:,2);params.limiter(1,2)]];
        end
        psi_sep = psi_wall;
        CQC = contourc(g.rg,g.zg,g.psirz-psiminimum,([1 1]*psi_sep-psiminimum)*0.999);
        fcc = find_closed_contour(CQC,params_fake);
        rbbbs = fcc.rb;
        zbbbs = fcc.zb;
        %ind = find(inpolygon(rbbbs,zbbbs,rwall,zwall)==1);
        %psi_wall = interp2(g.rg,g.zg,g.psirz,rwall(ind),zwall(ind),'cubic');
        %psi_sep = min(psi_wall);
        %CQC = contourc(g.rg,g.zg,g.psirz-psiminimum,([1 1]*psi_sep-psiminimum)*0.99);
        %fcc = find_closed_contour(CQC,params_fake);
        %rbbbs = fcc.rb;
        %zbbbs = fcc.zb;
    end
end

%% Check for snowflake divertor
if ~strcmp(plasma,'limited')
    if length(rS) > 1
        dxp = min(sqrt((rS(1) - rS(2:end)).^2 + (zS(1) - zS(2:end)).^2));
        if isfield(params,'geqdsk')
            s = dxp/params.geqdsk.amin;
        else
            s = dxp/params.a;
        end
        if s <= 0.3
            plasma = 'snowflake';
        end
    end
end

%% Normalizing psi
in = inpolygon(rE,zE,rbbbs,zbbbs);
fE = fE(in); rE = rE(in); zE = zE(in);
[fE,indfE] = max(fE);
rE = rE(indfE); zE = zE(indfE);
psi_maxis = fE;
g.psirzn = (g.psirz - psi_maxis)/(psi_sep - psi_maxis);
[rbmax,irbmax] = max(rbbbs);
zbmax = zbbbs(irbmax);
bpmid = interp2(g.rg,g.zg,sqrt(g.br.^2 + g.bz.^2),rbmax,zbmax,'spline');
lq = 0.63e-3/bpmid^1.19;
psisol = interp2(g.rg,g.zg,g.psirzn,linspace(rbmax,rbmax+lq,11),zbmax*ones(1,11),'cubic');
rbbbs = interp1(1:length(rbbbs),rbbbs,linspace(1,length(rbbbs),87));
zbbbs = interp1(1:length(zbbbs),zbbbs,linspace(1,length(zbbbs),87));

%% Sorting x-points by flux value
[fS,indS] = sort(fS);
rS = rS(indS); zS = zS(indS);

%% Output data
output.rg = g.rg;
output.zg = g.zg;
output.rbbbs = rbbbs;
output.zbbbs = zbbbs;
output.psirz = g.psirz;
output.psirzn = g.psirzn;
output.magflux = psi_maxis;
output.boundflux = psi_sep;
output.rE = rE;
output.zE = zE;
output.fE = fE;
output.rS = rS;
output.zS = zS;
output.fS = fS;
output.psisol = psisol;
output.plasma = plasma;

g.rmaxis = output.rE;
g.zmaxis = output.zE;
g.sibry = psi_sep;
g.simag = psi_maxis;
g.rbbbs = rbbbs;
g.zbbbs = zbbbs;
g.nbbbs = length(g.rbbbs);

writeg(g,'geqdsk_out')
g = readg('geqdsk_out');
g.filename = 'geqdsk';
g.time = 0;
g.shot = 0;
g.date     = datestr(now,'dd/mm/yyyy');
g.fittype = 'FBT';
delete('geqdsk_out')
writeg(g,'geqdsk_out')
output.geqdsk = readg('geqdsk_out');
delete('geqdsk_out')

end