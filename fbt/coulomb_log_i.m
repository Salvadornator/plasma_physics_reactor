function Ln = coulomb_log_i(ni,Ti,Zeff)

Ln = 30 - log(Zeff.^3.*sqrt(ni)./Ti.^1.5);

end