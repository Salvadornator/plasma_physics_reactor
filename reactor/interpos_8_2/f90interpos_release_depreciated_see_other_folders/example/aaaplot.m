% a.out > aaa.m
% >> aaa

set_defaults_matlab

yok=sin(out(:,1)).^2;
yokp=2.*sin(out(:,1)).*cos(out(:,1));
yokpp=2.*cos(out(:,1)).*cos(out(:,1)) - 2.*sin(out(:,1)).*sin(out(:,1));
yokint=out(:,1)./2 +-sin(2.*out(:,1))./4.;

taus=0.01;
[atau,btau,ctau,dtau]=interpos(43,in(:,1),in(:,2),out(:,1),taus,[0 0],[0 0]);

figure(1);clf
plotos(out(:,1),yok,'k-',[1 0]);
hold on
icol=0;
%icol=icol+1;plotos(out0(:,1),out0(:,2),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out(:,1),out(:,2),'-',0,0,colos(icol,:));
icol=icol+1;plotos(outper(:,1),outper(:,2),'--',0,0,colos(icol,:));
icol=icol+1;plotos(outsplibnd(:,1),outsplibnd(:,2),'--',0,0,colos(icol,:));
plotos(in(:,1),in(:,2),'k*',[0 6]);
%plot(out(:,1),at,'g-.')
legend('anal','spline','taus=0.01','per','splibnd 0.01')
ylabel('y')

figure(2);clf
plotos(out(:,1),yokp,'k-',[1 0]);
hold on
icol=0;
icol=icol+1;plotos(out0(:,1),out0(:,3),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out(:,1),out(:,3),'-',0,0,colos(icol,:));
icol=icol+1;plotos(outper(:,1),outper(:,3),'--',0,0,colos(icol,:));
icol=icol+1;plotos(outsplibnd(:,1),outsplibnd(:,3),'--',0,0,colos(icol,:));
legend('anal','spline','taus=0.01','per','splibnd 0.01')
ylabel('yprime')

figure(3);clf
plotos(out(:,1),yokpp,'k-',[1 0]);
hold on
icol=0;
icol=icol+1;plotos(out0(:,1),out0(:,4),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out(:,1),out(:,4),'-',0,0,colos(icol,:));
icol=icol+1;plotos(outper(:,1),outper(:,4),'--',0,0,colos(icol,:));
icol=icol+1;plotos(outsplibnd(:,1),outsplibnd(:,4),'--',0,0,colos(icol,:));
legend('anal','spline','taus=0.01','per','splibnd 0.01')
ylabel('yprimeprime')

figure(4);clf
plotos(out(:,1),yokint,'k-',[1 0]);
hold on
icol=0;
icol=icol+1;plotos(out0(:,1),out0(:,5),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out(:,1),out(:,5),'-',0,0,colos(icol,:));
icol=icol+1;plotos(outper(:,1),outper(:,5),'--',0,0,colos(icol,:));
icol=icol+1;plotos(outsplibnd(:,1),outsplibnd(:,5),'--',0,0,colos(icol,:));
legend('anal','spline','taus=0.01','per','splibnd 0.01')
ylabel('yint')
