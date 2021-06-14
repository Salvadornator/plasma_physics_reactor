% a.out > aaa.m
% >> aaa

set_defaults_matlab

figure(11);clf
errorbar(in(:,1),in(:,2),in(:,3),'k*');
hold on
icol=0;
icol=icol+1;plotos(out(:,1),out(:,2),'-',0,0,colos(icol,:));
ylabel('Te [eV]')
xlabel('\rho')
axis([0 1.2 0 6000])
legend('data','interpos')

figure(12);clf
errorbar(in(:,1),in(:,2),in(:,3),'k*');
hold on
icol=0;
icol=icol+1;plotos(out(:,1),out(:,2),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out0(:,1),out0(:,2),'-',0,0,colos(icol,:));
icol=icol+1;plotos(out02(:,1),out02(:,2),'--',0,0,colos(icol,:));
ylabel('Te [eV]')
xlabel('\rho')
axis([0 1.2 0 6000])

legend('data','interpos','std spline','std 2')

figure(13);clf
icol=0;
icol=icol+1;plotos(out(:,1),out(:,3),'-',0,0,colos(icol,:));
hold on
icol=icol+1;plotos(out0(:,1),out0(:,3),'-',0,0,colos(icol,:));
%icol=icol+1;plotos(out02(:,1),out02(:,3),'--',0,0,colos(icol,:));
ylabel('dTe/drho [eV]')
xlabel('\rho')
%axis([0 1.2 0 6000])

%legend('interpos','std spline','std 2')
legend('interpos','std spline')
