function struct = bootstrap_coefs(struct)

struct = nu_star_e(struct);
struct = nu_star_i(struct);
struct = trapped_fraction(struct);

% X = struct.ft_eff_33 and F33
X = struct.ft./(1 + (0.55-0.1*struct.ft).*sqrt(struct.nue) + 0.45*(1-struct.ft).*struct.nue./struct.Zeff.^1.5);
struct.bs_coeff.F33 = 1 - (1+0.36./struct.Zeff).*X + 0.59./struct.Zeff.*X.^2 - 0.23./struct.Zeff.*X.^3;

% struct.ft_eff_31 and L31
X = struct.ft./(1 + (1-0.1*struct.ft).*sqrt(struct.nue)+0.5*(1-struct.ft).*struct.nue./struct.Zeff);
struct.bs_coeff.L31 = (1 + 1.4./(struct.Zeff+1)).*X - 1.9./(struct.Zeff+1).*X.^2 + 0.3./(struct.Zeff+1).*X.^3;

% X = struct.ft_eff_32_ee and L32_ee | Y = struct.ft_eff_32_ei and L32_ei
X = struct.ft./(1 + 0.26*(1-struct.ft).*sqrt(struct.nue)+0.18*(1-0.37*struct.ft).*struct.nue./sqrt(struct.Zeff));
F32_ee =  (0.05+0.62*struct.Zeff).*(X-X.^4)./(struct.Zeff.*(1+0.44*struct.Zeff)) + (X.^2 - X.^4 - 1.20*(X.^3-X.^4))./(1+0.22*struct.Zeff) + 1.2./(1+0.5*struct.Zeff).*X.^4 - 1.2./(1+0.5*struct.Zeff);
Y = struct.ft./(1 + (1+0.6*struct.ft).*sqrt(struct.nue)+0.85*(1-0.37*struct.ft).*struct.nue.*(1+struct.Zeff));
F32_ei = -(0.56+1.93*struct.Zeff).*(Y-Y.^4)./(struct.Zeff.*(1+0.44*struct.Zeff)) + 4.95*(Y.^2 - Y.^4 - 0.55*(Y.^3-Y.^4))./(1+2.48*struct.Zeff) - 1.2./(1+0.5*struct.Zeff).*Y.^4 + 1.2./(1+0.5*struct.Zeff);
struct.bs_coeff.L32 = F32_ee + F32_ei;

% X = struct.ft_eff_34 and L34
X = struct.ft./(1 + (1-0.1*struct.ft).*sqrt(struct.nue)+0.5*(1-0.5*struct.ft).*struct.nue./struct.Zeff);
struct.bs_coeff.L34 = (1 + 1.4./(struct.Zeff+1)).*X - 1.9./(struct.Zeff+1).*X.^2 + 0.3./(struct.Zeff+1).*X.^3;

% alpha
a0 = -1.17*(1-struct.ft)./(1-0.22*struct.ft-0.19*struct.ft.^2);
struct.bs_coeff.alpha = ((a0+0.25*(1-struct.ft.^2).*sqrt(struct.nui))./(1+0.5*sqrt(struct.nui)) + 0.315*struct.nui.^2.*struct.ft.^6)./(1+0.15*struct.nui.^2.*struct.ft.^6);

end