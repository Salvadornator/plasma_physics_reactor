function M = mutual_inductance(r,z,G,R2,Z2,dR2,dZ2,N2,NR)

%% Defining filaments in coils 1 and 2
NZ2 = ceil(NR*dZ2/dR2);
if NZ2<1
    NZ2 = 1;
end
R2in = R2 - dR2/2 + dR2/NR/2;
R2out = R2 + dR2/2 - dR2/NR/2;
Z2in = Z2 - dZ2/2 + dZ2/NZ2/2;
Z2out = Z2 + dZ2/2 - dZ2/NZ2/2;
R2filaments = linspace(R2in,R2out,NR);
Z2filaments = linspace(Z2in,Z2out,NZ2);
[R2,Z2] = meshgrid(R2filaments,Z2filaments);

%% Calculating Green's functions between external coils and grid points
M = sum(sum(N2/NR/NZ2*interp2(r,z,G,R2,Z2,'spline')));

end