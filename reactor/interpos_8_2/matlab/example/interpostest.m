set_defaults_matlab

if ~exist('testcase')
  testcase='3';
end
testcase
switch testcase
  case '1'
    xx=[1 2 3 4 5 6];
    yy=[1.23 2.93 2.234 4.23 2.45 1.35 ];
    aa=interpos(13,xx,yy,[0:0.2:7]);
    figure(11)
    plot(xx,yy,'k*')
    hold on
    plot([0:0.2:7],aa,'r--')
    
    
    x=linspace(0,1,100);
    y=x.^2;
    [a,b]=interpos(13,x,y);

  case '2'

    a=[1:2000];b=a.^1.7;c=randn(1,14000000)*500000;tic;d=interpos(-43,a,b,c,3.5);toc
    %  (for timing but be careful plotting very large c, d)
    
  case '3'
    % periodic function with random noise
    x1=0; x2=2*pi;
    xlim1=-3;xlim2=8;
    x=linspace(x1,x2,30);
    xspl=linspace(-pi,x2+pi/2,1000);
    xspl=xspl(end:-1:1);
    epsilon=0.05;
    y=sin(x).^2+2.*(rand(1,length(x))-0.5).*epsilon;
    y_anal=sin(x).^2;
    yprime_anal=2.*sin(x).*cos(x);
    yprimeprime_anal=2.*(cos(x).*cos(x)-sin(x).^2);
    yint_anal=x./2-sin(2.*x)./4;
    y_analout=sin(xspl).^2;
    yprime_analout=2.*sin(xspl).*cos(xspl);
    yprimeprime_analout=2.*(cos(xspl).*cos(xspl)-sin(xspl).^2);
    yint_analout=xspl./2-sin(2.*xspl)./4;
    [yspl,yp,ypp,yint]=interpos(13,x,y,xspl);
    taus=7e-4;
    taus=2e-2;
    taus=0.01*mean(diff(x)).^3*100;
    taus=0.01
    % [atau,btau,ctau,dtau]=interpos(53,x,y,xspl,taus,[1 2],[0 0]);
    [atau,btau,ctau,dtau]=interpos(53,x,y,xspl,taus,[0 0],[0 0]);
    [atauper,btauper,ctauper,dtauper]=interpos(3,x(1:end),y(1:end),xspl,-1.,-1,x2-x1);
    figure(31);clf
    plotos(x,y,'k*',[0 6]);
    hold on
    %plotos(x,sin(x).^2,'k-',[1 0]);
    plotos(xspl,y_analout,'k-',[1 0]);
    icol=0;
    icol=icol+1;plotos(xspl,yspl,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,atau,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,atauper,'--',0,0,colos(icol,:));
    axis([xlim1 xlim2 -0.2 1.2])
    legend('data+eps','y anal','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)],2)
    title('y=sin(x+randn(1,length(x))*0.1).^2')
    
    figure(32);clf
    plotos(xspl,yprime_analout,'k-',[1 0]);
    hold on
    icol=0;
    icol=icol+1;plotos(xspl,yp,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,btau,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,btauper,'--',0,0,colos(icol,:));
    plotos(x,yprime_anal,'k*',[0 6]);
    title('yprime')
    legend('anal','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)],2)
    axis([xlim1 xlim2 -1.5 1.5])
    
    figure(33);clf
    plotos(xspl,yprimeprime_analout,'k-',[1 0]);
    hold on
    icol=0;
    icol=icol+1;plotos(xspl,ypp,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,ctau,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,ctauper,'--',0,0,colos(icol,:));
    plotos(x,yprimeprime_anal,'k*',[0 6]);
    title('yprimeprime')
    legend('anal','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)],2)
    axis([xlim1 xlim2 -5 5])
    
    figure(34);clf
    plotos(xspl,yint_analout,'k-',[1 0]);
    hold on
    icol=0;
    icol=icol+1;plotos(xspl,yint,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,dtau,'-',0,0,colos(icol,:));
    icol=icol+1;plotos(xspl,dtauper,'--',0,0,colos(icol,:));
    % plotos(x,yint_anal,'k*',[0 6]);
    title('yint')
    legend('anal','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)],2)
    axis([xlim1 xlim2 -3. 5])
    
  case '4'
    % take data as Te(rho)

    % $$$         load profdata_34175.mat
    % $$$         Tedata=profdata_plotted.meanprof_eff.te.data;
    % $$$         Testd=profdata_plotted.meanprof_eff.te.err;
    % $$$         Terho=profdata_plotted.meanprof_eff.te.rhoaxis_vol;
    % $$$         rhofit=profdata_plotted.fit.Xrhoaxis_vol;
    % $$$         ndata=length(Tedata);
    % $$$         nrhofit=length(rhofit);
    % save -ascii Teprof.data ndata Tedata Testd Terho nrhofit rhofit
    load Teprof.data; n1=Teprof(1);Tedata=Teprof(2:n1+1);Testd=Teprof(n1+2:2*n1+1);Terho=Teprof(2*n1+2:3*n1+1);
    n2=Teprof(3*n1+2);rhofit=Teprof(3*n1+3:3*n1+2+n2);
    figure(41);clf
    herr=errorbar(Terho,Tedata,Testd,'k*');
    hold on
    icol=0;
    % $$$         icol=icol+1;plotos(rhofit,profdata_plotted.fit.te,'-',[],[],colos(icol,:));
    taus=1e-6;
    [a,b,c,d]=interpos(13,Terho,Tedata,rhofit,taus,[],[],Testd);
    c(1)
    icol=icol+1;plotos(rhofit,a,'--',[],[],colos(icol,:));
    [a,b,c,d]=interpos(53,Terho,Tedata,rhofit,taus,[],[],Testd);
    c(1)
    icol=icol+1;plotos(rhofit,a,'--',[],[],colos(icol,:));
    [a,b,c,d]=interpos(3,[0; Terho],[Tedata(1); Tedata],rhofit,taus,[],[],[1.*Testd(1);Testd]);
    c(1)
    icol=icol+1;plotos(rhofit,a,'-',[],[],colos(icol,:));
    [a,b,c,d]=interpos(3,[0; Terho],[Tedata(1); Tedata],rhofit,taus,[0 0],[0 0],[1e4.*Testd(1);Testd]);
    c(1)
    icol=icol+1;plotos(rhofit,a,'-',[],[],colos(icol,:));
    [a,b,c,d]=interpos(3,[0; Terho],[Tedata(1); Tedata],rhofit,taus,[1 0],[0 0],[1e4.*Testd(1);Testd]);
    c(1)
    b(1)
    icol=icol+1;plotos(rhofit,a,'-',[],[],colos(icol,:));
    
    legend(num2str([1:icol+1]'))
    axis([0 1.2 0 6500])
    
  case '5'
    % magnetic data: time, amplitude, freq, approx_width
    cats=load('53290_catsOK.txt');
    figure(51);clf
    icol=0;
    icol=icol+1;plotos(cats(:,1),cats(:,2),'-',0,0,'k');
    icol=0;
    taus0=mean(diff(cats(:,1))).^3;
    [b0,db0]=interpos(13,cats(:,1),cats(:,2),taus0);
    hold on
    icol=icol+1;plotos(cats(:,1),b0,'--',0,0,'c');
    taus=taus0*10000;
    [b1,db1]=interpos(13,cats(:,1),cats(:,2),taus);
    icol=icol+1;plotos(cats(:,1),b1,'--',0,0,colos(icol,:));

    figure(52);clf
%    subplot(2,1,1);
    x =[59.75, 61.67, 62.48, 63.39, 65.97, 66.40];
    indxx=find(b1>0);
    indx(1)=indxx(1);
    for i=1:length(x)
      [zz ix]=min(abs(cats(:,1)-x(i)));
      indx(i+1)=ix;
    end
    plotos(cats(:,1),cats(:,2),'-',[1.5 0],0,'k');
    hold on
    for i=1:length(indx)-1
      plotos(cats(indx(i):indx(i+1),1),b1(indx(i):indx(i+1)),'-',0,0,colos(i,:));
      if i==1; hold on; end
    end
    figure(53);clf
%    subplot(2,1,2);
    for i=1:length(x)
      [zz ix]=min(abs(cats(:,1)-x(i)));
      indx(i+1)=ix;
    end
    for i=1:length(indx)-1
      %      plotos(sqrt(abs(b0)),db0./b0,'--',0,0,'c');
      plotos(sqrt(abs(b1(indx(i):indx(i+1)))),db1(indx(i):indx(i+1))./b1(indx(i):indx(i+1)),'-',0,0,colos(i,:));
      if i==1; hold on; end
    end
    axis([0.   0.03 -10.  15])
    
  case '6'
    % EQDSK
    fname='EQDSK.29866t1.5003';
    fname='eqdsksigns.36151t0.5000';
    fname='eqdsksigns.31837t1.0000';
    eqdskval=readeqdsk(fname);
    NR=3*length(eqdskval.rmesh);
    NZ=3*length(eqdskval.zmesh);
    rmeshfit=linspace(eqdskval.rmesh(1),eqdskval.rmesh(end),NR);
    zmeshfit=linspace(eqdskval.zmesh(1),eqdskval.zmesh(end),NR);
    tension1D=2e-6;
    tension1D=0.5.*((min(diff(eqdskval.rmesh))+min(diff(eqdskval.zmesh)))./2).^3
    ioptos=43;
    clear psioz_1 dpsidZ_1 d2psidZ2_1
    clear psioz0_1 dpsidZ0_1 d2psidZ20_1
    tic
    for iR=1:length(eqdskval.rmesh)
      [psiz_1(iR,:),dpsidZ_1(iR,:),d2psidZ2_1(iR,:)]= ...
          interpos(ioptos,eqdskval.zmesh,eqdskval.psi(iR,:),eqdskval.zmesh,tension1D);
      [psiz0_1(iR,:),dpsidZ0_1(iR,:),d2psidZ20_1(iR,:)]= ...
          interpos(ioptos,eqdskval.zmesh,eqdskval.psi(iR,:),eqdskval.zmesh);
    end
    toc
    tic
    clear psir_1 dpsidR_1 d2psidR2_1 dpsidZ_2 d2psidZdR_2 dpsidR_2 d2psidRZ2_2 
    clear psir0_1 dpsidR0_1 d2psidR20_1
    for iZ=1:length(eqdskval.zmesh)
      [psir_1(:,iZ),dpsidR_1(:,iZ),d2psidR2_1(:,iZ)]=interpos(ioptos,eqdskval.rmesh',eqdskval.psi(:,iZ),eqdskval.rmesh',tension1D);
      [dpsidZ_2(:,iZ),d2psidZdR_2(:,iZ)]=interpos(ioptos,eqdskval.rmesh',dpsidZ_1(:,iZ),eqdskval.rmesh',tension1D);
      [psir0_1(:,iZ),dpsidR0_1(:,iZ),d2psidR20_1(:,iZ)]=interpos(ioptos,eqdskval.rmesh',eqdskval.psi(:,iZ),eqdskval.rmesh');
    end
    for iR=1:length(eqdskval.rmesh)
      [dpsidR_2(iR,:),d2psidRZ2_2(iR,:)]=interpos(ioptos,eqdskval.zmesh,dpsidR_1(iR,:),eqdskval.zmesh,0);
    end
    toc
    % compute Grad-Shafranov terms
    % find pprime and ttprime on r,z mesh
    tic
    clear pprimeRZ FFprimeRZ
    for iR=1:length(eqdskval.rmesh)
      pprimeRZ(iR,:)= interpos(13,eqdskval.psimesh,eqdskval.pprime,(eqdskval.psi(iR,:)-eqdskval.psiaxis)./(eqdskval.psiedge-eqdskval.psiaxis));
      FFprimeRZ(iR,:)= interpos(13,eqdskval.psimesh,eqdskval.FFprime,(eqdskval.psi(iR,:)-eqdskval.psiaxis)./(eqdskval.psiedge-eqdskval.psiaxis));
    end
    toc
    ij=find((eqdskval.psi-eqdskval.psiaxis)./(eqdskval.psiedge-eqdskval.psiaxis)>1);
    pprimeRZ(ij)=0.;
    FFprimeRZ(ij)=0.;
    mu0=4e-7.*pi;
    RRR=eqdskval.rmesh'*ones(1,length(eqdskval.zmesh));
    rjphi=-RRR.^2.*mu0.*pprimeRZ-FFprimeRZ;
    gradshaf= d2psidR2_1 -1./RRR.*dpsidR_1 + d2psidZ2_1;
    gradshaf0= d2psidR20_1 -1./RRR.*dpsidR0_1 + d2psidZ20_1;
    % 2D spline 0 smoothing
    do2D=0
    if do2D
      sp0 = spaps({eqdskval.rmesh,eqdskval.zmesh},eqdskval.psi,0);  
      deriv0R=fnder(sp0,[1 0]);
      deriv0R2=fnder(sp0,[2 0]);
      deriv0Z2=fnder(sp0,[0 2]);
      [xx,yy] = ndgrid(eqdskval.rmesh,eqdskval.zmesh);
      dfdRall0=reshape(fnval(deriv0R,[xx(:) yy(:)]'),length(eqdskval.rmesh),length(eqdskval.zmesh));
      d2fdR2all0=reshape(fnval(deriv0R2,[xx(:) yy(:)]'),length(eqdskval.rmesh),length(eqdskval.zmesh));
      d2fdZ2all0=reshape(fnval(deriv0Z2,[xx(:) yy(:)]'),length(eqdskval.rmesh),length(eqdskval.zmesh));
      gradshaf2D0= d2fdR2all0 -1./RRR.*dfdRall0 + d2fdZ2all0;
    end
    [zzz ir]=min(abs(eqdskval.rmesh-eqdskval.raxis));
    [zzz iz]=min(abs(eqdskval.zmesh-eqdskval.zaxis));
    figure(61);clf
    plot(eqdskval.rmesh,rjphi(:,iz),'k')
    hold on
    plot(eqdskval.rmesh,gradshaf0(:,iz),'b--')
    plot(eqdskval.rmesh,gradshaf(:,iz),'r--')
    if do2D; plot(eqdskval.rmesh,gradshaf2D0(:,iz),'c:'); end
    xlabel(['R at Z=' num2str(eqdskval.zmesh(iz),'%.2f')])
    legend('R j_{\phi}','\Delta^* \psi std spline','\Delta^* \psi \tau spline')
    if do2D; legend('R j_{\phi}','\Delta^* \psi std spline','\Delta^* \psi \tau spline','2D 0spline'); end
    
    figure(62);clf
    plot(eqdskval.zmesh,rjphi(ir,:),'k')
    hold on
    plot(eqdskval.zmesh,gradshaf(ir,:),'r--')
    xlabel(['Z at R=' num2str(eqdskval.rmesh(ir),'%.2f')])
    legend('R j_{\phi}','\Delta^* \psi')
    
    figure(63);clf
    contour(eqdskval.rmesh,eqdskval.zmesh,eqdskval.psi',60,'k')
    hold on
    contour(eqdskval.rmesh,eqdskval.zmesh,psiz_1',60,'--')
    axis equal

    figure(64);clf
    plot(eqdskval.zmesh,d2psidZ20_1(ir,:),'-')
    hold on
    plot(eqdskval.zmesh,d2psidZ2_1(ir,:),'r-')
    plot(eqdskval.zmesh,d2psidR20_1(ir,:),'--')
    plot(eqdskval.zmesh,d2psidR2_1(ir,:),'r--')
    xlabel(['Z at R=' num2str(eqdskval.rmesh(ir),'%.2f')])
    legend('d^2\psi/dZ^2 std spline (\tau=0)','d^2\psi/dZ^2 \tau spline', ...
	   'd^2\psi/dR^2 std spline','d^2\psi/dR^2 \tau spline')
    
    figure(65);clf
    plot(eqdskval.rmesh,d2psidZ20_1(:,iz),'-')
    hold on
    plot(eqdskval.rmesh,d2psidZ2_1(:,iz),'r-')
    plot(eqdskval.rmesh,d2psidR20_1(:,iz),'--')
    plot(eqdskval.rmesh,d2psidR2_1(:,iz),'r--')
    xlabel(['R at Z=' num2str(eqdskval.zmesh(iz),'%.2f')])
    legend('d^2\psi/dZ^2 std spline (\tau=0)','d^2\psi/dZ^2 \tau spline', ...
	   'd^2\psi/dR^2 std spline','d^2\psi/dR^2 \tau spline')
    
  otherwise
  
    disp('not defined')
end
