function p = iterate_pedestal_reactor
%% Input parameters
pathname = '/home/fmsalvador/matlab/plasma_physics_project/reactor/scenarios';
shot     = 100;
eqtime   = 100;
Ip       = 8e6;
Paux     = 25e6;
Te0      = 11.5;
R0       = 2.8;
B0       = -6;

%% Calculating pedestal structure
Ti0      = Te0*1.15;
ne0      = 0.4*Ip/1e6/pi/0.9^2;
neped    = 0.75*ne0;
Zeff     = 1.5;
fz       = 5e-4; % Concentration of iron in the plasma
Teped    = Te0/3;   % Initial guess
Tiped    = Teped*0.85; % Initial guess
nesep    = 0.1*neped;
Tesep    = 0.1*Teped;
d_e      = 0.07;    % Initial guess
d_i      = d_e*1.2; % Initial guess

% Creating run directory
pathname = fullfile(pathname,sprintf('shot%06u',shot));
if ~exist(pathname,'dir')
    eval(['!mkdir ' pathname])
end
cd(pathname)
    
create_profiles_iterative(ne0,Te0,Ti0,nesep,Tesep,Teped,Tiped,neped,d_e,d_i,Zeff,fz,pathname);
p          = create_shape_reactor;
p.shot     = shot;
p.eqtime   = eqtime;
p.i_save   = 1;
p.jbs_fac  = 1;
p.d_e      = d_e;
p.d_i      = d_i;
p.bcp      = 0;
p.Ip       = Ip;
p.Bcentr   = B0;
p.Rcentr   = R0;
p.Paux     = Paux;
p.kinprofs = pathname;
p          = fbt(p);

fprintf('   --- STPED converged ---\n')

end