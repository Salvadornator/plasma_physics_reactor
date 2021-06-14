set_defaults_matlab

x=linspace(0,10,51);
x2=linspace(-1,11,61);
xspl=linspace(-1,11,601);
%xspl=linspace(4.1,4.3,11);
taus=1e-3;
y=(x-2).^1 .* (x-8).^2;
yp=1.*(x-2).^0 .* (x-8).^2 + (x-2).^1 .* 2.*(x-8).^1;
y2=(x2-2).^1 .* (x2-8).^2; % (x-2)*(x^2-16*x+64) = x^3-18*x^2+96*x-128
y2p=1.*(x2-2).^0 .* (x2-8).^2 + (x2-2).^1 .* 2.*(x2-8).^1;
ypp=0.*(x-2).^1 .* (x-8).^3 + 1.*(x-2).^0 .* 4.*(x-8).^1 + (x-2).^1 .* 2.*(x-8).^0;
y2pp=0.*(x2-2).^1 .* (x2-8).^3 + 1.*(x2-2).^0 .* 4.*(x2-8).^1 + (x2-2).^1 .* 2.*(x2-8).^0;
y2int=x2.^4./4-18*x2.^3/3+96.*x2.^2/2-128.*x2;
if ~exist('iopt')
  iopt=93
end
tic;[a,b,c, d]=interpos(iopt,x,y,xspl,taus,[0 0],[ypp(1) ypp(end)]);toc
%tic;[a,b,c, d]=interpos(x,y,xspl,60.*taus,[0 0],[-20 0],1+100.*(10.-x).^2);toc
%tic;[a,b,c, d]=interpos(x,y,xspl,60.*taus,[0 0],[-20 0],1+100.*(x).^2);toc
which interpos

figure(1);clf
plot(x,y,'k*-')
hold on
plot(x2,y2,'gs')
plot(xspl,a,'r--')
axis([-1 11 -300 100])

figure(2);clf
plot(x,yp,'k*-')
hold on
plot(x2,y2p,'gs-')
plot(xspl,b,'r--')
axis([-1 11 -50 150])

figure(3);clf
plot(x2,y2pp,'ks-')
hold on
plot(xspl,c,'r--')
axis([-1 11 -50 40])

figure(4);clf
plot(x2,y2int,'ks-')
hold on
plot(xspl,d,'r--')
axis([-1 11 -300 200])
