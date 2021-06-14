function [D, f, dr, dz]=create_D(r_vessel, z_vessel, psi_vessel, J_vessel)
%matrix
%    | D  -I  0 |
%A = |-I   D -I |
%    | 0  -I  D |


dr = r_vessel(2) - r_vessel(1);
dz = z_vessel(2) - z_vessel(1);
dz_dr2 = (dz / dr)^2;

% For central difference:
a = dz_dr2 * r_vessel ./ (r_vessel + (dr / 2));
b = dz_dr2 * r_vessel ./ (r_vessel - (dr / 2));
% For forward difference:
% a = dz_dr2 * (1 - (dr / r_vessel(2:end)))
% b = dz_dr2 * ones(r_vessel(2:end).size)
a = a(2:end);
b = b(2:end);
e = - 2 - a - b;


mu0=1.2566370614359173e-06;
% source in poisson equation
f = -mu0 * r_vessel * J_vessel * (dz^2);


% boundary conditions
f(1, :) = f(1, :) + psi_vessel(1, :);
f(length(r_vessel)-1, :) = f(length(r_vessel)-1, :) + psi_vessel(length(r_vessel), :);
f(:, 2) = f(:, 2) + psi_vessel(:, 1);
f(:, length(r_vessel)-1) = f(:, length(r_vessel)-1) + psi_vessel(:, length(r_vessel));

f = f(2:end, 2:end);

% first iteration D matrix
D = zeros(r.size - 2, r.size - 2);
[m,n]=size(D)

for i=1:m
    D(i, i) = - e(i);
    if 0 < i < m - 1
        D(i, i - 1) = - b(i);
        D(i, i + 1) = - a(i);
    elif i == 0
        D(i, i + 1) = - a(i);
    elif i == m - 1
        D(i, i - 1) = - b(i);
    end
end

