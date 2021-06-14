function G = green_em(Ro,Zo,R,Z,mode)
R = abs(R);
Mo = pi*4E-7;
h = Z - Zo;
u = sqrt((R + Ro).^2 + h.^2);
k = sqrt(4*Ro*R./u.^2);
d = sqrt((R - Ro).^2 + h.^2);
v = sqrt(R.^2 + Ro.^2 + h.^2);
w = sqrt(Ro.^2 - R.^2 - h.^2);
[K,E] = ellipke(k.^2);

switch mode
    case {'psi','mutual'}
        G = Mo*u.*((1 - 1/2*k.^2).*K - E);
    case 'br'
        G = Mo/2/pi*h./R./u.*(v.^2./d.^2.*E - K);
    case 'bz'
        G = Mo/2/pi./u.*(w.^2./d.^2.*E + K);
    case {'dpsidr','dMdr'}
        G = 2*pi*R.*(Mo/2/pi./u.*(w.^2./d.^2.*E + K));
    case {'dpsidz','dMdz'}
        G = -2*pi*R.*(Mo/2/pi*h./R./u.*(v.^2./d.^2.*E - K));
    case {'dbrdz','dbzdr'}
        G = Mo/2/pi./R./d.^2./u.^3.*((v.^2.*u.^2 - h.^2.*(d.^2 + k.^2.*u.^4./d.^2)).*E + (h.^2.*v.^2 - u.^2.*d.^2).*K);
    case 'dbzdz'
        G = Mo/2/pi*h./d.^2./u.^3.*(-(3*u.^2 + 4*v.^2.*w.^2./d.^2).*E + w.^2.*K);
    case 'dbrdr'
        G = -(Mo/2/pi*h./R./u.*(v.^2./d.^2.*E - K)./R + Mo/2/pi*h./d.^2./u.^3.*(-(3*u.^2 + 4*v.^2.*w.^2./d.^2).*E + w.^2.*K));
end
G(R == Ro & Z == Zo) = 0;

end









