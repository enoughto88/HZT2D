%Initial condition of phi and te

%clear all
close all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [5 5 5.5 4.5]);


PI = 3.1415926535893;
nx = 240;
ny = 240;
nxd= 24;
nyd= 24;
XL = 2.0;
dx = XL/(nx-1);
xx = 0:dx:XL;
ratio = 10;
tecen1  = zeros(1,nx);
tecen2  = zeros(1,nx);
tecen3  = zeros(1,nx);
tecen4  = zeros(1,nx);
tecen5  = zeros(1,nx);
tecen6  = zeros(1,nx);
tecen7  = zeros(1,nx);
tedata = zeros(nxd*nyd,1);

pos = 1.0;
Tep = 20.0;
Tem = 4.0;
sig = 0.20;
for i = 1:1:nx
   tecen1(1,i) = (Tep-Tem)*exp(-1*(xx(i)-pos)^2/(2*sig^2))+Tem;
   tecen2(1,i) = (Tep-Tem)*exp(-1*(xx(i)-pos+0.1)^2/(2*sig^2))+Tem;
   tecen3(1,i) = (Tep-Tem)*exp(-1*(xx(i)-pos+0.2)^2/(2*sig^2))+Tem;
   tecen4(1,i) = (Tep-Tem)*exp(-1*(xx(i)-pos+0.3)^2/(2*sig^2))+Tem;
   tecen5(1,i) = (Tep+5-Tem)*exp(-1*(xx(i)-pos)^2/(2*sig^2))+Tem;
   tecen6(1,i) = (Tep+10-Tem)*exp(-1*(xx(i)-pos)^2/(2*sig^2))+Tem;
   tecen7(1,i) = (Tep+15-Tem)*exp(-1*(xx(i)-pos)^2/(2*sig^2))+Tem;
end


for j = 1:1:nyd
   for i = 1:1:nxd
      ip = (i-1)*ratio+1;
      jp = (j-1)*ratio+1;
      tedata(nxd*(j-1)+i) = tecen1(1,ip);
   end
end
%MFdata  = horzcat(phdata,tedata);
%dlmwrite(strcat('/Users/kawashima/Dropbox/HZT2D/program/','initialphite.dat'),MFdata, 'delimiter', '\t','precision','%20.10e');


%plot(xx,tecen1,'k-',xx,tecen2,'b-',xx,tecen3,'r-',xx,tecen4,'g-','LineWidth',1.5);
plot(xx*25,tecen4,'k-',xx*25,tecen1,'b--','LineWidth',1.5);
grid off
set(gca, 'XLim', [0,XL*25]);
set(gca, 'YLim', [0,25]);
set(gca,'XTick',[0.0 10 20.0 30 40 50])
set(gca,'XTickLabel',{'0.0','10.0','20.0','30.0','40.0','50.0'})
set(gca,'YTick',[0 5 10 15 20 25 30 35 40])
set(gca,'YTickLabel',{'0.0','5.0','10.0','15.0','20.0','25.0','30.0','35.0','40.0'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Electron temperature, eV','FontSize',16)
%legend('p = 1.0','p = 0.9','p = 0.8','p = 0.7','Location','NorthEast')
legend('Original','Shifted','Location','NorthEast')
saveas(figure(1),'/Users/kawashima/Dropbox/HZT2D/figure/datafile/gaussTe1.png');


