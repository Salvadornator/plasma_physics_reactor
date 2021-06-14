function out = readg(filename,varargin)
% This MATLAB routine reads GEQDSK files. It can also be used to plot the
% data in the GEQDSK file.
%
%   Example: g = readg(Filename[,machine])
%
%            Filename          : A string with the full path to the GEQDSK file
%            Machine (optional): It can be 'nstx', 'nstxu' or 'd3d'
%
%   g = readg('/home/canal/efit/nstx/132543/g132543.00719_641')
%
%   g = readg('/p/tsc/gpaganin/projects/ncc_design/free_boundary/2.0MA/equilibrium/geqdsk','nstxu')
%
%   g = readg('/home/canal/efit/nstx/132543/shot132543_kinetic_NB/g132543.00719_641','nstx')
%
%   g = readg('/home/canal/matlab/equil_reconst/scenarios/shot000700/g000700.00100','tcabr');

out.filename = filename;
if ~(exist(filename,'file') == 2)
    out = [];
    disp('   --- Warning: File does not exist ---')
    return
end

if nargin == 1
    iplot = 0;
else
    iplot = 1;
    machine = varargin{1};
end

%% Reading GEQDSK file
fid = fopen(filename);
    dum = textscan(fid,'%s',4);
    shot = dum{1}{3};
    time = dum{1}{4};
    out.shot = str2double(shot(2:end));
    out.time = str2double(time(1:end-2))/1000;
    out.fittype = dum{1}{1};
    out.date = dum{1}{2};
    dum = textscan(fid,'%d',1);
    if double(dum{1}) >= 999
        dum = num2str(dum{1});
        out.idum = str2double(dum(1));
        out.nw = str2double(dum(2:5));
        out.nh = str2double(dum(6:9));
    else
        out.idum = double(dum{1});
        dum = textscan(fid,'%d',1);
        out.nw = double(dum{1});
        dum = textscan(fid,'%d',1);
        out.nh = double(dum{1});
    end
    dum = textscan(fid,'%f',5);
    out.rdim = double(dum{1}(1));
    out.zdim = double(dum{1}(2));
    out.rcentr = double(dum{1}(3));
    out.rleft = double(dum{1}(4));
    out.zmid = double(dum{1}(5));
    out.rg = linspace(out.rleft,out.rleft+out.rdim,out.nw).';
    out.zg = linspace(-out.zdim/2 - out.zmid,out.zdim/2 + out.zmid,out.nw).';
    dum = textscan(fid,'%f',5);
    out.rmaxis = double(dum{1}(1));
    out.zmaxis = double(dum{1}(2));
    out.simag = double(dum{1}(3));
    out.sibry = double(dum{1}(4));
    out.bcentr = double(dum{1}(5));
    dum = textscan(fid,'%f',10);
    out.current = double(dum{1}(1));
    dum = textscan(fid,'%f',out.nw);
    out.fpol = double(dum{1});
    out.Ipol = 2*pi/4/pi/1e-7*abs(out.fpol - out.fpol(1));
    %out.rcentr = out.fpol(end)/out.bcentr;
    dum = textscan(fid,'%f',out.nw);
    out.pres = double(dum{1});
    dum = textscan(fid,'%f',out.nw);
    out.ffprim = double(dum{1});
    dum = textscan(fid,'%f',out.nw);
    out.pprime = double(dum{1});
    dum = textscan(fid,'%f',out.nw.*out.nh);
    out.psirz = reshape(double(dum{1}),out.nw,out.nh).'; % Stream function: minus the poloidal flux per radian ( -psi_pol / 2 pi )
    dum = textscan(fid,'%f',out.nw);
    out.R = [];
    out.Z = [];
    out.psirzn = (out.psirz - out.simag)/(out.sibry - out.simag); % Normalized stream function
    out.thetarz = [];
    out.jacobian = [];
    out.br = [];
    out.bz = [];
    out.dbrdr = [];
    out.dbrdz = [];
    out.dbzdr = [];
    out.dbzdz = [];
    out.jr = [];
    out.jz = [];
    out.jphi = [];
    out.qpsi = double(dum{1});
    dum = textscan(fid,'%d',2);
    out.nbbbs = double(dum{1}(1));
    out.limitr = double(dum{1}(2));
    dum = textscan(fid,'%f',2*out.nbbbs);
    dum = reshape(double(dum{1}),2,out.nbbbs);
    out.rbbbs = double(dum(1,:)).';
    out.zbbbs = double(dum(2,:)).';
    dum = textscan(fid,'%f',2*out.limitr);
    dum = reshape(double(dum{1}),2,out.limitr);
    out.rlim = double(dum(1,:)).';
    out.zlim = double(dum(2,:)).';
    dum = textscan(fid,'%d',1);
    out.kvtor = double(dum{1});
    if isempty(out.kvtor)
        out.kvtor = 0;
    end
    dum = textscan(fid,'%f',1);
    out.rvtor = double(dum{1});
    if isempty(out.rvtor)
        out.rvtor = 1.70000005;
    end
    dum = textscan(fid,'%d',1);
    out.nmass = double(dum{1});
    if isempty(out.nmass)
        out.nmass = 0;
    end
    if out.kvtor > 0
        dum = textscan(fid,'%f',out.nw);
        out.pressw = double(dum{1});
        dum = textscan(fid,'%f',out.nw);
        out.pwprim = double(dum{1});
    end
    if out.nmass > 0
        dum = textscan(fid,'%f',out.nw);
        out.dmion = double(dum{1});
    end
    dum = textscan(fid,'%f',out.nw);
    if ~isempty(dum{1})
        out.rhovn = double(dum{1});    
    end
fclose(fid);
out.psin = linspace(0,1,out.nw).';
out.amin = (max(out.rbbbs) - min(out.rbbbs))/2;

% Calculating flux surface averages
out = flux_average(out);
out.nG = out.current/1e6/pi/out.amin^2*1e20;

%% Plotting GEQDSK file data
if iplot
    [rbmax,ind_rbmax] = max(out.rbbbs);
    zbmax = out.zbbbs(ind_rbmax);
    bpmid = interp2(out.rg,out.zg,sqrt(out.br.^2 + out.bz.^2),rbmax,zbmax,'spline');
    lq = 0.63e-3/bpmid^1.19;
    psi_sol = interp2(out.rg,out.zg,out.psirzn,linspace(rbmax,rbmax+lq,11),zbmax*ones(1,11),'cubic');
    c = mycolor;
    % Loading machine structure
    switch machine
        case 'd3d'
            tok_data_struct = load('/home/canal/matlab/machines/d3d_obj_mks_struct_129129.mat');
        case 'nstx'
            tok_data_struct = load('/home/canal/matlab/machines/nstx_obj_12nov_6565.mat');
        case 'nstxu'
            tok_data_struct = load('/home/canal/matlab/machines/nstxu_obj_15feb_129129.mat');
        case 'east'
            tok_data_struct = load('/home/canal/matlab/machines/east_obj_2016_3333.mat');
        case ['tca','tcabr']
            tok_data_struct = [];
            options.iblackbg = 0;
    end
    
    if ~strcmp(machine,'tca') && ~strcmp(machine,'tcabr')
        tok_data_struct = tok_data_struct.tok_data_struct;
    end

    % Plotting
    nfig = figure(287);
    clf
    set(nfig,'Name','Tokamak Equilibrium')
    subplot(1,3,1);
    contour(out.rg,out.zg,out.psirzn,0.1:0.1:0.9,'color',c(2,:))
    hold on
    contour(out.rg,out.zg,out.psirzn,linspace(1,5,51),'color',c(2,:))
    contour(out.rg,out.zg,out.psirzn,psi_sol,'color',c(2,:))
    contour(out.rg,out.zg,out.psirzn,[1 1],'color',c(3,:),'linewidth',2)
    plot(out.rbbbs,out.zbbbs,'color',c(3,:),'linewidth',2)
    plot(out.rmaxis,out.zmaxis,'+k','linewidth',2,'markersize',10)
    if strcmp(machine,'tcabr') && ~strcmp(machine,'tca')
        plot_tcabr('fig',nfig)
    elseif strcmp(machine,'tca')
        plot_tcabr('fig',nfig,'limiter',1)
    else
        plot_tok_geo(tok_data_struct)
    end
    hold off
    axis equal
    xlabel('R (m)')
    ylabel('Z (m)')
    grid on
    title({['#' num2str(out.shot) ' @ ' num2str(out.time*1000) ' ms'];['\lambda_q = ' sprintf('%3.1f',lq*1000) ' mm']})
    if strcmp(machine,'nstx') || strcmp(machine,'nstxu')
        axis([0 2.2 -2 2])
    elseif strcmp(machine,'d3d')
        axis([0.5 3.2 -2.1 2.1])
    elseif strcmp(machine,'tcabr') || strcmp(machine,'tca')
        axis([0.25 1.0 -0.7 0.7])
    end
    
    h(1) = subplot(2,3,2);
    plot(out.psin,out.pprime/1e6*out.rcentr,'color',c(3,:),'linewidth',3)
    hold on
    plot(out.psin,out.ffprim/(pi*4e-7)/1e6/out.rcentr,'color',c(1,:),'linewidth',3)
    xlabel('\Psi_N')
    grid on
    legend('R_0 P '' ( MA / m^2 )','FF '' / \mu_0 R_0 ( MA / m^2 )','location','northwest')
    legend('boxoff')
    %axis([0 1 0 Inf])

    h(2) = subplot(2,3,3);
    plot(out.psin,out.pres/1e3,'color',c(3,:),'linewidth',3)
    xlabel('\Psi_N')
    ylabel('Pressure ( kPa )')
    grid on
    xlim([0 Inf])

    h(3) = subplot(2,3,5);
    plot(out.psin,out.qpsi,'color',c(3,:),'linewidth',3)
    hold on
    plot(out.psin,out.s,'color',c(1,:),'linewidth',3)
    plot(out.psin,out.alpha,'color',c(6,:),'linewidth',3)
    hold off
    grid on
    xlabel('\Psi_N')
    legend('q','s','\alpha','location','northwest')
    legend('boxoff')
    axis([0 1 0 6])
    
    h(4) = subplot(2,3,6);
    plot(out.psin,out.jpar/1e6,'color',c(3,:),'linewidth',3)
    hold on
    plot(out.psin,out.jtor/1e6,'color',c(1,:),'linewidth',3)
    plot(out.psin,out.jpol/1e6,'k','linewidth',3)
    plot(out.psin,out.Ipol/1e6,'color',c(6,:),'linewidth',3)
    hold off
    grid on
    %ylim([0 5])
    xlabel('\Psi_N')
    legend('J_{||} ( MA / m^2 )','J_{tor} ( MA / m^2 )','J_{pol} ( MA / m^2 )','I_{pol} ( MA )','location','northeast')
    legend('boxoff')
    
    linkaxes(h,'x')
end

%% Function 1
function outt = get_contour_levels(CCC,xo,yo)
Npoints = 0;
kk = 0;
iii = 0;
while kk+Npoints < length(CCC)
    kk = kk + Npoints + 1;
    Npoints = CCC(2,kk);
    if CCC(1,kk+1) == CCC(1,kk+Npoints) && CCC(2,kk+1) == CCC(2,kk+Npoints)
        [in,~] = inpolygon(xo,yo,CCC(1,kk+1:kk+Npoints),CCC(2,kk+1:kk+Npoints));
        if in == 1
            iii = iii + 1;
            Rb{iii} = CCC(1,kk+1:kk+Npoints);
            Zb{iii} = CCC(2,kk++1:kk+Npoints);
            V(iii) = CCC(1,kk);
            N(iii) = CCC(2,kk);
        end
    end
end
outt.C = CCC;
if exist('V','var')
    outt.value = V;
else
    outt = [];
    fprintf('\n       --- Warning: No closed contour was found ---\n')
    return
end
outt.n_points = N;
outt.Rb = Rb;
outt.Zb = Zb;
end

%% Function 2
function f = dFdx(x,F)
if length(x) == length(F)
    Lx = length(x);
    f = zeros(size(x));
    for iii = 2:Lx-1
        %f(iii) = (-F(iii+2) + 8*F(iii+1) - 8*F(iii-1) + F(iii-2))/(x(iii+2) - x(iii-2))/3;
        f(iii) = (F(iii+1) - F(iii-1))/(x(iii+1) - x(iii-1));
    end
    f([1 Lx]) = interp1(x(2:Lx-1),f(2:Lx-1),x([1 Lx]),'pchip');
else
    error('Vectors must have same dimension')
end
end

%% Function 3
function col = mycolor
% Defining colors
    col = zeros(8,3);
    col(1,:) = [237 32 36];   % Red
    col(2,:) = [164 194 219]; % Light blue
    col(3,:) = [57 83 164];   % Dark Blue
    col(4,:) = [50 180 80];   % Light Green
    col(5,:) = [0 0 0];       % Black
    col(6,:) = [0 110 0];     % Green
    col(7,:) = [256 0 256];   % Magenta
    col(8,:) = [0 256 256];   % Magenta
    col = col/256;
end

%disp('Done!')

end