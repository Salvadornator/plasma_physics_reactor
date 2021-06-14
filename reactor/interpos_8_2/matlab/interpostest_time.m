x=linspace(0,6*pi,50);
x=linspace(0,6*pi,50000);
xspl2=linspace(0.1,5*pi,3000000);
xspl2=linspace(-3*pi,9*pi,30);
xspl2=linspace(-3*pi,9*pi,3000000);
xspl3=linspace(0.1,5.9*pi,3000000);
taus=1e-6;
tic;[atauper,btauper,ctauper,d1]=interpos(13,x(1:end-1),sin(x(1:end-1)).^2,xspl2,taus,-1,6.*pi);toc
tic;[atauper,btauper,ctauper,d]=interpos(13,x(1:end),sin(x(1:end)).^2,xspl2,taus,[2 2],[0 0]);toc
tic;[atauper,btauper,ctauper,d1]=interpos(3,x(1:end-1),sin(x(1:end-1)).^2,xspl3,taus,-1,6.*pi);toc
tic;[atauper,btauper,ctauper,d]=interpos(3,x(1:end),sin(x(1:end)).^2,xspl3,taus,[2 2],[0 0]);toc
which interpos
