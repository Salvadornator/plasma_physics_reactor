function params = solve_GS(params,iplot)
    
  % Solving Grad-Shafranov equation
  fprintf('       --- Solving Grad-Shafranov equation\n')
    [RJcurr,ZJcurr]      = meshgrid(params.jcur.rg,params.jcur.zg);
    struct_temp.r        = params.rmesh.';
    struct_temp.z        = params.zmesh.';
    struct_temp.psi      = -2*pi*(params.psi - params.psic - params.psioh - params.psiv);
    %struct_temp.psicoils = params.psic + params.psioh + params.psiv;
    if ~isempty(params.jplasma)
        struct_temp.J = params.jplasma;
    else
        [Rcurr,Zcurr]   = meshgrid(params.rmesh,params.zmesh);
        struct_temp.J = interp2(params.jcur.rg,params.jcur.zg,params.jcur.jplasma,Rcurr,Zcurr);
    end
    try
        save('./run_equil.mat','-struct','struct_temp')
    catch ME
        if (strcmp(ME.identifier,'MATLAB:save:couldNotWriteFile'))
            msg = 'FBT is trying to save a file but you are running from a directory that is not yours. Change directory to be able to run the code.';
            ME = MException('MATLAB:save:couldNotWriteFile',msg);
       end
       throw(ME)
    end
    system('/opt/anaconda5/bin/python3.6 /home/canal/matlab/equil_reconst/fbt/gs_solver/grad_shafranov.py');
    delete('./run_equil.mat')
    fields = load('./rec_fields.mat');
    delete('./BCR_*.npy')
    if exist('./fort.0','file') == 2
        delete('./fort.0')
    end
    if exist('./fort.6','file') == 2
        delete('./fort.6')
    end
    delete('./rec_fields.mat')
    params.jcur.jplasma = interp2(params.rmesh,params.zmesh,fields.J_rec,RJcurr,ZJcurr);
    if isfield(params,'geqdsk')
        in = inpolygon(RJcurr,ZJcurr,params.geqdsk.rbbbs,params.geqdsk.zbbbs);
        params.jcur.jplasma(~in) = 0;
        params.jcur.jplasma = params.jcur.jplasma/sum(params.jcur.jplasma(:)*mean(diff(params.jcur.rg))*mean(diff(params.jcur.zg)))*params.Ip;
    end

    params.psi = -fields.psi_rec/2/pi + params.psic + params.psioh + params.psiv;
    fields     = find_xpoints(params);
    if isempty(fields)
        params = [];
        return
    end
    params.psin         = fields.psirzn;    % Normalized poloidal flux distribution
    params.fbt.rE       = fields.rE;      % R of the Elliptic points
    params.fbt.zE       = fields.zE;      % Z of the Elliptic points
    params.fbt.fE       = fields.fE;      % Poloidal flux at the elliptic points
    params.fbt.rS       = fields.rS;      % R of the Saddle points
    params.fbt.zS       = fields.zS;      % Z of the Saddle points
    params.fbt.fS       = fields.fS;      % Poloidal flux at the saddle points
    params.Rbb          = fields.rbbbs;
    params.Zbb          = fields.zbbbs;
    params.geqdsk       = fields.geqdsk;
    params.psisol       = fields.psisol;

    % Plotting total poloidal flux distribution
  if iplot
    c = mycolor;
    if strcmp(fields.plasma,'limited')
        icol = [3 1 6 12 7];
    elseif strcmp(fields.plasma,'snowflake')
        icol = [1 3 6 12 7];
    else % Single-null or double-null
        icol = [3 1 6 12 7];
    end

    figure(1)
    clf
    psinnn = params.psin;
    [rnnn,znnn] = meshgrid(params.rmesh,params.zmesh);
    innn = inpolygon(rnnn,znnn,params.geqdsk.rbbbs,params.geqdsk.zbbbs);
    psinnn(~innn) = NaN;
    contour(params.rmesh,params.zmesh,psinnn,0.1:0.1:0.9,'color',c(2,:))
    hold on
    if strcmp(params.gridtype,'large')
        contour(params.rmesh,params.zmesh,params.psin,linspace(1,max(params.psin(:)),161),'color',c(2,:))
    else
        contour(params.rmesh,params.zmesh,params.psin,linspace(1,8,81),'color',c(2,:))
    end
    contour(params.rmesh,params.zmesh,params.psin,fields.psisol,'color',c(2,:))
    plot(params.Rb,params.Zb,'ko','color',[1 1 1]*0.7,'linewidth',2)
    plot(params.Rbfix,params.Zbfix,'ko','color',[1 1 1]*0.7,'linewidth',2,'markersize',10)

    if strcmp(fields.plasma,'limited')
        plot(params.Rbb,params.Zbb,'color',c(3,:),'linewidth',2)
    end
    if strcmp(fields.plasma,'snowflake')
        vec = length(params.fbt.fS):-1:1;
    else
        vec = 1:length(params.fbt.fS);
    end
    if ~strcmp(fields.plasma,'limited')
        for ixp = 1:length(vec)
            contour(params.rmesh,params.zmesh,params.psi,[1 1]*params.fbt.fS(vec(ixp)),'color',c(icol(ixp),:),'linewidth',2)
        end
    end
    if params.Rcentr == 0.62
    	plot_tcabr
    else
        plot_reactor
    end
    plot(params.fbt.rE,params.fbt.zE,'+k','linewidth',2,'markersize',10)
    plot(params.Rmaxis,params.Zmaxis,'+','color',[1 1 1]*0.7,'linewidth',1,'markersize',10)
    if ~strcmp(fields.plasma,'limited')
        plot(params.Rxp,params.Zxp,'x','color',[1 1 1]*0.7,'linewidth',2,'markersize',12)
        plot(params.fbt.rS,params.fbt.zS,'xk','linewidth',2,'markersize',12)
    end
    plot(params.Rsp,params.Zsp,'^k','linewidth',2,'markersize',6)
    if params.i_vessel
        colorbar
    end
    
    % Plotting vessel currents
    if params.i_vessel && params.iter > 1 && std(params.jv(:,end)) ~= 0
        Ivmax = max(abs(params.jv(:,end)/1e3));
        col = b2r_colormap(-Ivmax,Ivmax);
        colormap(col)
        dIv = 2*Ivmax/(size(col,1)-1);
        for iv = 1:size(params.jv,1)
            icol = 1 + floor((params.jv(iv,end)/1e3 + Ivmax)/dIv);
            fill([params.vv(iv,1)-params.dRv(iv)/2 params.vv(iv,1)+params.dRv(iv)/2 params.vv(iv,1)+params.dRv(iv)/2 params.vv(iv,1)-params.dRv(iv)/2 params.vv(iv,1)-params.dRv(iv)/2],[params.vv(iv,2)-params.dZv(iv)/2 params.vv(iv,2)-params.dZv(iv)/2 params.vv(iv,2)+params.dZv(iv)/2 params.vv(iv,2)+params.dZv(iv)/2 params.vv(iv,2)-params.dZv(iv)/2],col(icol,:),'edgecolor','none')
        end
        plot(params.W(:,1),params.W(:,2),'k')
        h = colorbar;
        ylabel(h,'J_{s,vv} ( kA / m )')
    end
    %axis([0 0.9 -0.3 0.3])

    Ithreshold = 6;
    figure(2)
    clf
    bar(params.Ia/1e3,'k','BarWidth',0.5)
    hold on
    if params.Rcentr == 0.62
        Ccritic = find(abs(params.Ia)/1e3>Ithreshold);
    else
        Ccritic = find(abs(params.Ia)/1e3>Ithreshold*1e6);
    end
    Icrit = zeros(1,length(params.Ia)); Icrit(Ccritic) = params.Ia(Ccritic);
    if ~isempty(Ccritic)
        bar(Icrit/1e3,'r','BarWidth',0.5)
    end
    if params.Rcentr == 0.62
        plot([0 17],[1 1]*Ithreshold,'k--')
        plot([0 17],-[1 1]*Ithreshold,'k--')
    end
    hold off
    %axis([0 17 -8 8])
    ylabel('I_{PF} ( kA )')
    xlabel('PF Coil #')
    drawnow
    
  end

end