function struct = bootstrap_j(struct)

struct = bootstrap_coefs(struct);
gp  = struct.bs_coeff.L31./struct.pe.*dFdx(struct.psin,struct.ptot)/(struct.sibry-struct.simag);
gte = struct.bs_coeff.L32./struct.Te.*dFdx(struct.psin,struct.Te)/(struct.sibry-struct.simag);
gti = struct.bs_coeff.L34.*struct.bs_coeff.alpha.*(struct.ptot - struct.pe)./struct.pe./struct.Ti.*dFdx(struct.psin,struct.Ti)/(struct.sibry-struct.simag);
struct.j_bs = -struct.rcentr*struct.pe.*(gp + gte + gti);
struct.j_bs_gp = -struct.rcentr*struct.pe.*gp;
struct.j_bs_gte = -struct.rcentr*struct.pe.*gte;
struct.j_bs_gti = -struct.rcentr*struct.pe.*gti;

end