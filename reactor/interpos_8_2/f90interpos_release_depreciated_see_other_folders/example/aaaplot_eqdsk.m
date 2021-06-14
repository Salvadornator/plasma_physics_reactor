% a.out > aaa.m
% >> aaa

set_defaults_matlab
set(0,'DefaultaxesFontName','a14');
set(0,'DefaulttextFontName','a14');
set(0,'DefaultUicontrolFontName','a14')
set(0,'FixedWidthFontName','a14')
set(0,'DefaultaxesFontSize',10);
set(0,'DefaulttextFontSize',10);
set(0,'DefaultaxesXGrid','off');
set(0,'DefaultaxesYGrid','off');


nrin=length(rmesh);
nzin=length(zmesh);
nrout=length(rmeshout);
nzout=length(zmeshout);

figure(161);clf
psieff=NaN*ones(nrin,nzin);
for iz=1:nzin
  psieff(:,iz)=psi((iz-1)*nrin+[1:nrin])';
end
[haa,hbb]=contour(rmesh,zmesh,psieff',100,'k-');
hlevel=get(hbb,'LevelList');
hold on
psiouteff=NaN*ones(nrout,nzout);
for iz=1:nzout
  psiouteff(:,iz)=psiout((iz-1)*nrout+[1:nrout])';
end
contour(rmeshout,zmeshout,psiouteff',hlevel,'--')
legend('orig','spline with tension')
ylabel('\psi(R,Z) contours')
axis equal

figure(162);clf
irmid=fix(nrout/2);
izmid=fix(nzout/2);
rjphieff=NaN*ones(nrout,nzout);
for iz=1:nzout
  rjphieff(:,iz)=rjphi((iz-1)*nrout+[1:nrout])';
end
plot(rmeshout,rjphieff(:,izmid),'k-')
hold on
gradshafeff=NaN*ones(nrout,nzout);
for iz=1:nzout
  gradshafeff(:,iz)=gradshaf((iz-1)*nrout+[1:nrout])';
end
plot(rmeshout,gradshafeff(:,izmid),'-')
legend('R j_{phi}','\Delta*')
ylabel(['values at z=' num2str(zmeshout(izmid),'%.2f')])

figure(163);clf
plot(rmeshout,rjphieff(irmid,:),'k-')
hold on
plot(rmeshout,gradshafeff(irmid,:),'-')
legend('R j_{phi}','\Delta*')
ylabel(['values at R=' num2str(rmeshout(irmid),'%.2f')])

figure(164);clf
set(gcf,'pos',[962   254   512 551])
set(gcf,'paperpositionmode','auto')
plot(rzplas(:,1),rzplas(:,2),'k-')
hold on
plot(rzplasout0(:,1),rzplasout0(:,2),'b-')
plot(rzplasout(:,1),rzplasout(:,2),'r-')
legend('orig','std spline','with tension',4)
ylabel('plasma boundary')
axis equal
grid off

figure(165);clf
plot(thetarhoin(:,1),thetarhoin(:,2),'k-')
hold on
plot(thetarhopppout0(:,1),thetarhopppout0(:,2),'b-')
plot(thetarhopppout(:,1),thetarhopppout(:,2),'r-')
legend('orig','std spline','with tension')
ylabel('plasma boundary rho(theta)')
xlabel('\theta')

figure(166);clf
plot(thetarhopppout0(:,1),thetarhopppout0(:,3),'b-')
hold on
plot(thetarhopppout(:,1),thetarhopppout(:,3),'r-')
legend('std spline','with tension')
ylabel('d\rho / d\theta')
xlabel('\theta')

figure(167);clf
plot(thetarhopppout0(:,1),thetarhopppout0(:,4),'b-')
hold on
plot(thetarhopppout(:,1),thetarhopppout(:,4),'r-')
legend('std spline','with tension')
ylabel('d^2\rho / d\theta^2')
xlabel('\theta')

figure(168);clf
set(gcf,'pos',[962   254   512 551])
set(gcf,'paperpositionmode','auto')
subplot(3,1,1)
plot(thetarhoin(:,1),thetarhoin(:,2),'k-')
hold on
plot(thetarhopppout0(:,1),thetarhopppout0(:,2),'b-')
plot(thetarhopppout(:,1),thetarhopppout(:,2),'r-')
hl1=legend('orig','std spline','with tension',0);
ylabel('\rho(\theta)')
set(hl1,'pos',[0.38   0.8  0.20 0.10])

subplot(3,1,2)
plot(thetarhopppout0(:,1),thetarhopppout0(:,3),'b-')
hold on
plot(thetarhopppout(:,1),thetarhopppout(:,3),'r-')
legend('std spline','with tension',2)
ylabel('d\rho / d\theta')

subplot(3,1,3)
plot(thetarhopppout0(:,1),thetarhopppout0(:,4),'b-')
hold on
plot(thetarhopppout(:,1),thetarhopppout(:,4),'r-')
legend('std spline','with tension',3)
ylabel('d^2\rho / d\theta^2')
xlabel('\theta')

modsubplot([0 2*pi],'grid off',0)
