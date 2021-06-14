% a.out > aaa.m
% >> aaa

set_defaults_matlab

taus=0.01;
[a0 b0 c0]=interpos(13,in(:,1),in(:,2),out(:,1),0.);
[at bt ct]=interpos(13,in(:,1),in(:,2),out(:,1),taus);
yok=sin(out(:,1)).^2;
yokp=2.*sin(out(:,1)).*cos(out(:,1));
yokpp=2.*cos(out(:,1)).*cos(out(:,1)) - 2.*sin(out(:,1)).*sin(out(:,1));
yokint=out(:,1)./2 + sin(2.*out(:,1))./4.;

figure(1);clf
plot(in(:,1),in(:,2),'k*')
hold on
plot(out(:,1),yok,'k-')
plot(out(:,1),out(:,2),'r--')
plot(outsplibnd(:,1),outsplibnd(:,2),'m-')
plot(out(:,1),at,'g-.')

figure(2);clf
plot(out(:,1),yokp,'k-')
hold on
plot(out(:,1),b0,'b')
plot(out(:,1),out(:,3),'r-')
plot(outsplibnd(1:60:end,1),outsplibnd(1:60:end,3),'ro')
plot(out(:,1),bt,'g-.')

figure(3);clf
plot(out(:,1),yokpp,'k-')
hold on
plot(out(:,1),c0,'b')
plot(out(:,1),out(:,4),'r-')
plot(outsplibnd(1:60:end,1),outsplibnd(1:60:end,4),'ro')
plot(out(:,1),ct,'g-.')
