function struct = j_ohm(struct)

spt = spitzer(struct.ne,struct.Te,struct.Zeff);
struct.sigma_spt_par = 1./spt.eta_par;
struct.sigma_spt_perp = 1./spt.eta_perp;
struct.sigma_neo_par = struct.sigma_spt_par.*struct.bs_coeff.F33;
struct.sigma0 = struct.sigma_neo_par.*struct.R02_R2;
struct.I_bs = struct.fpol(end)/struct.rcentr.*trapz(struct.psin,struct.Vp./struct.fpol.*struct.j_bs*(struct.sibry-struct.simag));
%struct.I_ps = 1/struct.rcentr*trapz(struct.psin,struct.Vp.*struct.j_ps.*(struct.sibry-struct.simag));
%struct.Vl = 2*pi*struct.rcentr^2*(struct.current - struct.jbs_fac*struct.I_bs - struct.jbs_fac*struct.I_ps)/struct.fpol(end)/trapz(struct.psin,struct.Vp./struct.fpol.*struct.sigma0.*(struct.sibry-struct.simag));
%struct.j_ohm = struct.sigma0.*struct.Vl/2/pi/struct.rcentr;
struct.I_sigma0 = struct.fpol(end)/2/pi/struct.rcentr^2*trapz(struct.psin,struct.Vp./struct.fpol.*struct.sigma0*(struct.sibry-struct.simag));

end