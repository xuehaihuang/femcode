clear;

h=1:8;
h=(1/2).^h;

eta=[1.0245E-01	3.0067E-02	1.1924E-02	6.2114E-03	3.2139E-03	1.7375E-03	1.0300E-03	6.5592E-04];
zeta=[6.4552E-02	1.6856E-02	5.7318E-03	2.9844E-03	1.7388E-03	9.9585E-04	5.7149E-04	3.3753E-04];
etatilde=[1.0231E-01	2.9848E-02	1.1748E-02	6.1283E-03	3.1755E-03	1.7165E-03	1.0165E-03	6.4661E-04];
zetatilde=[6.4323E-02	1.6461E-02	5.3556E-03	2.8073E-03	1.6669E-03	9.5878E-04	5.4673E-04	3.1905E-04];


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


line([1 1 2 1]+3.5,[-4 -5 -5 -4]-3.3,'LineWidth',2,'Color',[0 0 0]);

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

