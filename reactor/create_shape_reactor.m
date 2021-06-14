function params = create_shape_reactor

%% Plasma shape parameters
params          = get_fbt_params_reactor;
params.Nb       = 100;
params.Rmaxis   = 2.7;
params.Zmaxis   = 0.0;
params.Rcentr   = 2.8;
params.a        = 0.9;
params.kappa    = 2.1;
params.deltau   = 0.3;
params.deltal   = 0.9;
params.sqrur    = 0; 
params.sqrul    = 0;
params.sqrll    = 0;
params.sqrlr    = 0;
params.Rxp      = 1.85;
params.Zxp      = -1.8;
params.plasmaid = 'diverted';
params.imode    = 'hmode';
params.kinprofs = '/home/fmsalvador/matlab/plasma_physics_project/reactor/scenarios';

params = get_boundary_points(params);

% Plotting boundary
figure(1)
clf
plot_reactor
plot(params.Rb,params.Zb,'k.','linewidth',2,'markersize',8)
hold off
title('Plasma Boundary')
drawnow

end