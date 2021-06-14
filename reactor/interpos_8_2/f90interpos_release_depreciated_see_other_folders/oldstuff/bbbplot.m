set_defaults_matlab

% bbb

taus=0.01;
[a0 b0 c0]=interpos(13,in(:,1),in(:,2),out(:,1),0.);
[at bt ct]=interpos(13,in(:,1),in(:,2),out(:,1),taus);
yok=sin(out(:,1)).^2;
yokp=2.*sin(out(:,1)).*cos(out(:,1));
yokpp=2.*cos(out(:,1)).*cos(out(:,1)) - 2.*sin(out(:,1)).*sin(out(:,1));

figure(1);clf
plot(in(:,1),in(:,2),'k*')
hold on
plot(out(:,1),yok,'k-')
plot(out(:,1),out(:,2),'r--')
plot(outsplibnd(1:60:end,1),outsplibnd(1:60:end,2),'ro')
plot(out(:,1),at,'g-.')
plot(outper(:,1),outper(:,2),'c--')
plot(outper(:,1)+6*pi,outper(:,2),'c--')
plot(outper(:,1)-6*pi,outper(:,2),'c--')

figure(2);clf
plot(out(:,1),yokp,'k-')
hold on
plot(out(:,1),b0,'b')
plot(out(:,1),out(:,3),'r-')
plot(outsplibnd(1:60:end,1),outsplibnd(1:60:end,3),'ro')
plot(out(:,1),bt,'g-.')
plot(outper(:,1),outper(:,3),'c--')
plot(outper(:,1)+6*pi,outper(:,3),'c--')
plot(outper(:,1)-6*pi,outper(:,3),'c--')

figure(3);clf
plot(out(:,1),yokpp,'k-')
hold on
plot(out(:,1),c0,'b')
plot(out(:,1),out(:,4),'r-')
plot(outsplibnd(1:60:end,1),outsplibnd(1:60:end,4),'ro')
plot(out(:,1),ct,'g-.')
plot(outper(:,1),outper(:,4),'c--')
plot(outper(:,1)+6*pi,outper(:,4),'c--')
plot(outper(:,1)-6*pi,outper(:,4),'c--')

figure
plot(out(:,1),yokpp,'k-')
hold on
plot(outper(:,1),outper(:,4),'r--')
