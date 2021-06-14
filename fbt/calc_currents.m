function params = calc_currents(params)

% --- Lagrange solver - P, A and b Matrices
fprintf('       --- Calculating coil currents\n')
rho_cu = 1.68e-8;
R = rho_cu*params.coils(:,5)*2*pi.*params.coils(:,1)./(params.coils(:,3).*params.coils(:,4)./params.coils(:,5));
R = diag([R; 0]);
Wp = params.wp*diag([ones(1,params.Nc_E) 0.25*ones(1,params.Nc_F+params.Nc_D) 0]);
kk = 0;
DD = zeros(size(params.coils,1)^2-size(params.coils,1),size(params.coils,1)+1);
for ii = 1:size(params.coils,1)
    for jj = 1:size(params.coils,1)
        if ii ~= jj
            kk = kk + 1;
            if mod(jj,17)==0
                DD(kk,[ii jj]) = [0 0];
            else
                DD(kk,[ii jj]) = [1 -1];%/((params.coils(ii,1) - params.coils(jj,1))^2 + (params.coils(ii,2) - params.coils(jj,2))^2);
            end
        end
    end
end
P = 3e-10*Wp.*R + 2e-15*params.wd*(DD.'*DD)/2;
A = [-2*pi*params.fbt.Gbc, -ones(length([params.Rb params.Rxp params.Rsf params.Rbfix]),1)];
b = -(params.fbt.psibp + params.fbt.Gboh*params.Ioh(end) + params.fbt.Gbv*params.Iv(:,end))*(-2*pi);
Wb = eye(length([params.Rb params.Rxp params.Rsf params.Rbfix]));

%% Physical connections between coils
%Mconnections = zeros(1,17);
%Mconnections(1,[15 16]) = [1 -1];
%Mconnections(2,[15 16]) = [1  1];
%Mconnections(1,[5 9]) = [1 1];
%Mconnections(1,[5 9]) = [1 1];
%Mconnections(2,4:5) = [1 -1];
%Mconnections(2,[4 12]) = [1 1];
%Mconnections(7,7:8) = [1 -1];
%Mconnections(5,9:10) = [1 -1];
%Mconnections(6,11:12) = [1 -1];
Mconnections = [];

%% Volt-second constrain
if params.i_voltsec
    Mvoltsec = -2*pi*[params.fbt.Gcmaxis, 0];
else
    Mvoltsec = [];
end

%% Constraint equations
C = [[params.fbt.Gbfixc; params.fbt.Gxpc; params.fbt.Gsfc; params.fbt.Gspc]*(-2*pi), -ones(size([params.fbt.Gbfixc; params.fbt.Gxpc; params.fbt.Gsfc; params.fbt.Gspc;],1),1)];
C = [C; [params.fbt.brcxp, zeros(size(params.fbt.brcxp,1),1)]];
C = [C; [params.fbt.bzcxp, zeros(size(params.fbt.bzcxp,1),1)]];
C = [C; [params.fbt.brcsf, zeros(size(params.fbt.brcsf,1),1)]];
C = [C; [params.fbt.bzcsf, zeros(size(params.fbt.bzcsf,1),1)]];
C = [C; [params.fbt.dbrdrcsf, zeros(size(params.fbt.dbrdrcsf,1),1)]];
%C = [C; [dbrdzcsf, zeros(size(dbrdzcsf,1),1)]];
C = [C; [params.fbt.dbzdrcsf, zeros(size(params.fbt.dbzdrcsf,1),1)]];
C = [C; [params.fbt.dbzdzcsf, zeros(size(params.fbt.dbzdzcsf,1),1)]];
C = [C; Mvoltsec];
C = [C; Mconnections];

d = -[[params.fbt.psibfixp+params.fbt.Gbfixoh*params.Ioh(end)+params.fbt.Gbfixv*params.Iv(:,end); params.fbt.psixpp + params.fbt.Gxpoh*params.Ioh(end) + params.fbt.Gxpv*params.Iv(:,end); params.fbt.psisfp + params.fbt.Gsfoh*params.Ioh(end) + params.fbt.Gsfv*params.Iv(:,end); params.fbt.psispp + params.fbt.Gspoh*params.Ioh(end) + params.fbt.Gspv*params.Iv(:,end)]*(-2*pi); ...
       params.fbt.brxpp + params.fbt.brxpoh*params.Ioh(end) + params.fbt.brxpv*params.Iv(:,end); ...
       params.fbt.bzxpp + params.fbt.bzxpoh*params.Ioh(end) + params.fbt.bzxpv*params.Iv(:,end); ...
       params.fbt.brsfp; ...
       params.fbt.bzsfp; ...
       params.fbt.dbrdrsfp; ...
       params.fbt.dbzdrsfp; ...
       params.fbt.dbzdzsfp; ...
      -params.volt_sec*ones(size(Mvoltsec,1),1); ...
       zeros(size(Mconnections,1),1)];
%d = -[psibfixp; psixpp; psisfp; psispp; brpx; bzpx; brpsf; bzpsf; dbrdrpsf; dbrdzpsf; dbzdrpsf; dbzdzpsf; zeros(size(Mconnections,1),1)];

%% Solve for coil currents
T = [[P + A.'*Wb*A; C], [C.'; zeros(size(C,1))]];
V = [A.'*Wb*b; d];
X = T\V;
if params.iter < 2
    params.Ia = X(1:params.Nc_E+params.Nc_F+params.Nc_D);
else
    params.Ia = params.iter_fac*params.Ia + (1 - params.iter_fac)*X(1:params.Nc_E+params.Nc_F+params.Nc_D);
end

%% Computing updated external coil poloidal flux distribution
fprintf('       --- Calculating poloidal flux due to coil currents\n')
params.psic = zeros(length(params.greens.r));
for icoil = 1:params.Nc_E
    eval(['params.psic = params.psic - params.greens.E' num2str(icoil) '.G/2/pi*params.Ia(' num2str(icoil) ');'])
end
for icoil = 1:params.Nc_F
    eval(['params.psic = params.psic - params.greens.F' num2str(icoil) '.G/2/pi*params.Ia(' num2str(icoil+params.Nc_E) ');'])
end
for icoil = 1:params.Nc_D
    eval(['params.psic = params.psic - params.greens.D' num2str(icoil) '.G/2/pi*params.Ia(' num2str(icoil+params.Nc_E+params.Nc_F) ');'])
end
fprintf('       --- Calculating poloidal flux due to ohmic current\n')
params.psioh = -params.greens.OH.G/2/pi*params.Ioh(end);
fprintf('       --- Calculating poloidal flux due to vessel currents\n')
[Rm,Zm] = meshgrid(params.rmesh,params.zmesh);
params.psiv = zeros(size(Rm));
for ibp = 1:size(params.vv,1)
    params.psiv = params.psiv - green_em(params.vv(ibp,1),params.vv(ibp,2),Rm,Zm,'mutual')*params.Iv(ibp,end)/2/pi;
    if 0
        figure(1);clf;plot_tcabr;
        contour(params.rmesh,params.zmesh,temp*params.Iv(ibp,end),101)
        colorbar; caxis([-1 1]*5e-6); drawnow
    end
end
fprintf('       --- Calculating total poloidal flux distribution\n')
params.psi = params.psip + params.psic + params.psioh + params.psiv;

end