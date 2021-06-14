function kinetic_profiles
% This function calculates the plasma kinetc profiles based on:
%  ne, Te, p_total and an assumed flat Zeff profile
%
%  Necessário usar MATLAB versão: 8.5.0.197613 (R2015a) para
%  conseguir visualizar os símbolos em greek font
%
%  Também é necessário "load path" para o meu diretório onde
%  está a função readg.m: /home/canal/matlab/diagnostics/efit/
%  Talvez, seja mais apenas rodar minha função startup para ajustar o
%  ambiente: /home/canal/matlab/startup.m
%
%

%% Assumed Zeff
Zeff = 2.5;

%% Loading kinetic profiles
g = readg('~/efit/tcabr/35047/newgeqdsk',0,'');
ne = load('~/efit/tcabr/35047/profile_ne.smoothed');
te = load('~/efit/tcabr/35047/profile_te.smoothed');
w = load('~/efit/tcabr/35047/profile_omega.smoothed');

%% Defining profiles on the same grid
psin = ne(:,1);
ne = ne(:,2)*1e20;
te = te(:,2)*1e3;
w = w(:,2)*1e3;
ptot = interp1(g.psin,g.pres,psin,'pchip');

%% Computing profiles
pe = 1.602e-19*ne.*te;
nc = ne*(Zeff - 1)/30;
ni = ne - 6*nc;
ti = (ptot - pe)./(ni + nc)/1.602e-19;
pion = 1.602e-19*ni.*ti;
pc = 1.602e-19*nc.*ti;

%% Plotting
c = mycolormap;
figure(27)
clf
plot(psin,ne/1e19,'color',c(1,:),'linewidth',3)
hold on
plot(psin,ni/1e19,'color',c(3,:),'linewidth',3)
plot(psin,6*nc/1e19,'color',c(5,:),'linewidth',3)
xlabel('\psi_N')
title('Number Density ( 1x10^{19} m^{-3} )')
legend('n_e','n_i','6 n_C')

figure(28)
clf
plot(psin,te,'color',c(1,:),'linewidth',3)
hold on
plot(psin,ti,'color',c(3,:),'linewidth',3)
xlabel('\psi_N')
title('Temperature ( eV )')
legend('T_e','T_i')

figure(29)
clf
plot(psin,w/1e3,'color',c(5,:),'linewidth',3)
xlabel('\psi_N')
title('Plasma Toroidal Rotation ( krad/s )')

figure(30)
clf
plot(psin,ptot/1e3,'color',c(2,:),'linewidth',3)
hold on
plot(psin,pe/1e3,'color',c(1,:),'linewidth',3)
plot(psin,pion/1e3,'color',c(3,:),'linewidth',3)
plot(psin,pc/1e3,'color',c(5,:),'linewidth',3)
xlabel('\psi_N')
title('Pressure ( kPa )')
legend('p_{total}','p_e','p_i','p_C')

%% Extra function
function c = mycolormap
    c = zeros(12,3);
    c(1,:) = [237 32 36];   % Red
    c(2,:) = [164 194 219]; % Light blue
    c(3,:) = [57 83 164];   % Dark Blue
    c(4,:) = [50 180 80];   % Light Green
    c(5,:) = [0 0 0];       % Black
    c(6,:) = [0 110 0];     % Green
    c(7,:) = [256 0 256];   % Magenta
    c(8,:) = [0 256 256];   % Cyan
    c(9,:) = [0 0 256];     % Normal blue
    c(10,:) = [256 0 0];    % Normal red 
    c(12,:) = [256 102 154];% Lilac
    c = c/256;
end

end