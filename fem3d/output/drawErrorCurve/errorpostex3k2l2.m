clear;

h=1:8;
h=(1/2).^h;

eta=[1.4581E-01	7.6200E-02	5.4883E-02	3.2316E-02	1.8222E-02	1.1103E-02	7.1767E-03	4.7797E-03];
zeta=[7.2751E-02	3.1236E-02	2.2980E-02	1.2469E-02	6.8334E-03	4.2848E-03	2.8642E-03	1.9498E-03];
etatilde=[1.4579E-01	7.6159E-02	5.4863E-02	3.2307E-02	1.8218E-02	1.1100E-02	7.1751E-03	4.7786E-03];
zetatilde=[7.2698E-02	3.1136E-02	2.2933E-02	1.2445E-02	6.8214E-03	4.2781E-03	2.8601E-03	1.9471E-03];


h=-log(h);
eta=log(eta);
zeta=log(zeta);
etatilde=log(etatilde);
zetatilde=log(zetatilde);

plot(h,eta,'o-',h,zeta,'<-',h,etatilde,'*-',h,zetatilde,'p-','LineWidth',2);
leg=legend('$\ln\eta_h(\mathbf{\sigma}_h, u_h, f)$','$\ln\zeta_h(\mathbf{\sigma}_h, u_h, f)$','$\ln\tilde{\eta}_h(\mathbf{\sigma}_h, u_h, f)$','$\ln\tilde{\zeta}_h(\mathbf{\sigma}_h, u_h, f)$');
%leg=legend('$\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$|u-u_h|_{1}$','$|\!|\!|u-u_h|\!|\!|$','Location','southwest');
set(leg,'Interpreter','latex');
set(xlabel('$\ln(1/h)$'),'Interpreter','latex');
%set(ylabel('$\ln$Errors'),'Interpreter','latex');

%set (gcf,'Position',[200, 100, 720, 540])
%set (gcf,'Position',[200, 100, 800, 600], 'color','w')
% reset: set (gcf,'Position',[232, 246, 560, 420], 'color','w')


line([1 1 2 1]+2.2,[-4 -5 -5 -4]-1,'LineWidth',2,'Color',[0 0 0]);

%line([1 1 2 1]+2.9,[-7 -9 -9 -7]-6.2,'LineWidth',2,'Color',[0 0 0]);
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

