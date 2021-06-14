
testcase='3'

switch testcase
      case '1'
        xx=[1 2 3 4 5 6];
        yy=[1.23 2.93 2.234 4.23 2.45 1.35 ];
        aa=interpos_f90all(13,xx,yy,[0:0.2:7]);
        figure
        plot(xx,yy,'k*')
        hold on
        plot([0:0.2:7],aa,'r--')


        x=linspace(0,1,100);
        y=x.^2;
        [a,b]=interpos_f90all(13,x,y);

      case '2'

        a=[1:2000];b=a.^1.7;c=randn(1,14000000)*500000;tic;d=interpos_f90all(-43,a,b,c,3.5);toc
        %  (for timing but be careful plotting very large c, d)

      case '3'

        x=linspace(0,6.*pi,100);
        xspl=linspace(0,6.*pi,1000);
        epsilon=0.1;
        y=sin(x+rand(1,length(x))*epsilon).^2;
        y=sin(x).^2+(rand(1,length(x))-0.5).*2.*epsilon;
        yprime=2.*sin(x).*cos(x);
        yprimeprime=2.*(cos(x).*cos(x)-sin(x).^2);
        [yspl,yp,ypp]=interpos_f90all(13,x,y,xspl);
        taus=7e-4;
        taus=2e-2;
        [atau,btau,ctau]=interpos_f90all(13,x,y,xspl,taus);
        [atauper,btauper,ctauper]=interpos_f90all(13,x(1:end-1),y(1:end-1),xspl,taus,-1,6.*pi);
        figure
        plot(x,y,'k*')
        hold on
        plot(x,sin(x).^2,'k-')
        plot(xspl,yspl)
        plot(xspl,atau,'r-')
        plot(xspl,atauper,'c--')
        legend('data','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)])
        title('y=sin(x+randn(1,length(x))*0.1).^2')

        figure
        plot(x,yprime,'k*')
        hold on
        plot(xspl,yp)
        plot(xspl,btau,'r-')
        plot(xspl,btau,'c--')
        title('yprime')
        legend('data','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)])

        figure
        plot(x,yprimeprime,'k*')
        hold on
        plot(xspl,ypp)
        plot(xspl,ctau,'r--')
        plot(xspl,ctauper,'c--')
        title('yprimeprime')
        legend('data','spline',['taus=' num2str(taus)],['per taus=' num2str(taus)])

end
