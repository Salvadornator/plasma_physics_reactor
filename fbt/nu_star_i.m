function struct = nu_star_i(struct)

struct.Lni = coulomb_log_i(struct.ni,struct.Ti,struct.Zeff);
struct.nui = 4.9e-18*struct.qpsi.*struct.rcentr.*struct.ni.*struct.Zeff.^4.*struct.Lni./(struct.Ti.^2.*(struct.amin./struct.rcentr).^1.5);

end