clear;

h=1:8;
h=(1/2).^h;

errorsigma0=[2.9659E-02	1.9320E-02	1.3769E-02	8.4860E-03	4.6893E-03	2.4479E-03	1.2469E-03	6.2866E-04];
erroru1=[3.7178E-03	1.5094E-03	5.7463E-04	1.8725E-04	5.3456E-05	1.4180E-05	3.6397E-06	9.2230E-07];
erroruE=[6.0563E-02	3.8481E-02	2.1237E-02	1.0578E-02	5.0971E-03	2.4836E-03	1.2263E-03	6.0974E-04];


h=-log(h);
errorsigma0=log(errorsigma0);
erroru1=log(erroru1);
erroruE=log(erroruE);

plot(h,errorsigma0,'o-',h,erroru1,'<-',h,erroruE,'*-','LineWidth',2);
leg=legend('$\ln\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$\ln|u-u_h|_{1}$','$\ln|\!|\!|u-u_h|\!|\!|$');
%leg=legend('$\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$|u-u_h|_{1}$','$|\!|\!|u-u_h|\!|\!|$','Location','southwest');
set(leg,'Interpreter','latex');
set(xlabel('$\ln(1/h)$'),'Interpreter','latex');
%set(ylabel('$\ln$Errors'),'Interpreter','latex');

%set (gcf,'Position',[200, 100, 720, 540])
%set (gcf,'Position',[200, 100, 800, 600], 'color','w')
% reset: set (gcf,'Position',[232, 246, 560, 420], 'color','w')


line([1 1 2 1]+3,[-4 -5 -5 -4]-2.5,'LineWidth',2,'Color',[0 0 0]);

line([1 1 2 1]+2.5,[-7 -9 -9 -7]-4,'LineWidth',2,'Color',[0 0 0]);
%line([3.5 2.5 3.5 3.5],[-3 -3 -4 -3],'LineWidth',2,'Color',[0 0 0]);

%line([2.5 2.5 3.5 2.5],[-1 -2 -2 -1],'LineWidth',2,'Color',[0 0 0]);

%line([6 6 8 6],[0 -1 -1 0],'LineWidth',2,'Color',[0 0 0]);



%plot(h,errorEk3l0,'p-',h,errorEk3l1,'<-',h,errorEk3l2,'o-',h,errorEk3l3,'s-','LineWidth',2);
%line([1.5 1.5 2.5 1.5],[-6 -8 -8 -6],'LineWidth',2,'Color',[0 0 0]);
%legend('l=0', 'l=1', 'l=2', 'l=3');

%subplot(2,1,1)
%plot(h,error1k2l1,'p-',h,error1k2l0,'<-',h,error1k1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');
%grid on;
%subplot(2,1,2)
%plot(h,errorEk2l1,'p-',h,errorEk2l0,'<-',h,errorEk1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');

grid on;

