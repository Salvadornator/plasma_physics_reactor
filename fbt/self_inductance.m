function L = self_inductance(R,S,N)

Muo = pi*4e-7;
a = sqrt(S/pi);
k = 2*sqrt(R.*abs(R-a))./(2*R-a);
[K,E] = ellipke(k.^2);
L = Muo*N.^2.*((2*R-a).*((1 - 1/2*k.^2).*K - E) + R/4);
%L = Muo*N.^2.*R.*(log(8*R./a) - 7/4);

end