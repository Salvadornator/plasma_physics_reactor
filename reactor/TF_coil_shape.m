function [R,Z,sigma,b] = TF_coil_shape(R1,R2,R0,B0)

% Calculating TF coil geometry                 
mu0 = 4e-7*pi;    % mu0 constant
B1 = B0*R0/R1;          % Magnetic Field Intensity at R=R1
p1 = B1^2/2/mu0;      % Pressure due to the Magnetic Field at R=R1
sigma = p1*R1*log(R2/R1); %Tension due Magnetic Field
b = sqrt(R2*R1);
k = 0.4999999*log(R1/R2);    % k squared

R = linspace(R1,R2,1001);
f = log(b./R)./sqrt(k^2-log(b./R).^2);
R = flip(R); f = flip(f);
Z = cumtrapz(R,f);

R = [flip(R) R R(end)];
Z = real([-flip(Z) Z -Z(end)]);

end