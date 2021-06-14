function struct = nu_star_e(struct)

struct.Lne = coulomb_log_e(struct.ne,struct.Te);
struct.nue = 6.921e-18*struct.qpsi.*struct.rcentr.*struct.ne.*struct.Zeff.*struct.Lne./(struct.Te.^2.*(struct.amin./struct.rcentr).^1.5);

end