function out = spitzer(ne,Te,Zeff)

e = 1.6022e-19; me = 9.1094e-31; eo = 8.8542e-12;
Fz = (1 + 1.198*Zeff + 0.222*Zeff.^2)./(1 + 2.966*Zeff + 0.753*Zeff.^2);
Ln = coulomb_log_e(ne,Te);
out.eta_par = 4*sqrt(2*pi)/3*e^2*sqrt(me)/(4*pi*eo)^2/e^1.5.*Zeff.*Fz.*Ln./Te.^1.5;
out.eta_perp = 1.96*out.eta_par;

end