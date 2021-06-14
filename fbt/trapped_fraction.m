function struct = trapped_fraction(struct)

% Calculating magnetic fields
dr = mean(diff(struct.rg));
dz = mean(diff(struct.zg));
[R,~] = meshgrid(struct.rg,struct.zg);
[dpsidr,dpsidz] = gradient(struct.psirz,dr,dz);
BR =  sign(struct.current)./R.*dpsidz;
BZ = -sign(struct.current)./R.*dpsidr;
BPOL = sqrt(BR.^2 + BZ.^2);
struct.br = BR; struct.bz = BZ;

% Getting contour lines
psin = [struct.psin(3:end-1); 0.999];
C = contourc(struct.rg,struct.zg,struct.psirzn,psin);
C_out = get_contour_levels(C,struct.rmaxis,struct.zmaxis);
if isempty(C_out)
    out = [];
    return
end

struct.fc = zeros(size(psin));
struct.ft = zeros(size(psin));
xx = linspace(0,1,129);
for ii = 1:length(C_out.n_points)
    dl = [0 sqrt((C_out.Rb{ii}(2:end) - C_out.Rb{ii}(1:end-1)).^2 + (C_out.Zb{ii}(2:end) - C_out.Zb{ii}(1:end-1)).^2)];
    l  = cumsum(dl);

    % Interpolating profiles
    bpol = interp2(struct.rg,struct.zg,BPOL,C_out.Rb{ii},C_out.Zb{ii},'cubic');
    bphi = interp1(linspace(0,1,struct.nh),struct.fpol,C_out.value(ii),'pchip')./C_out.Rb{ii};

    % Fraction of circulating and trapped particles
    B2 = trapz(l,(bpol.^2 + bphi.^2)./bpol)/trapz(l,1./bpol);
    Bs = sqrt(bpol.^2 + bphi.^2);
    Bmax = max(sqrt(bpol.^2 + bphi.^2));
    msqrt = zeros(size(xx));
    for jj = 1:length(xx)
        msqrt(jj) = trapz(l,sqrt(1 - xx(jj)*Bs/Bmax)./bpol)/trapz(l,1./bpol);
    end
    struct.fc(ii) = 3/4*B2/Bmax^2*trapz(xx,xx./msqrt);
    struct.ft(ii) = 1 - struct.fc(ii);
end
struct.fc  = interp1(C_out.value,struct.fc,struct.psin,'pchip');
struct.ft  = interp1(C_out.value,struct.ft,struct.psin,'pchip');

end