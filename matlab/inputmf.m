%Magnetic field calculater of SPT-100
%
%clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');

PI = 3.1415926535893;
nx = 240;
ny = 240;
nxd= 48;
nyd= 48;
ratio = 5;
XL = 0.05;
YL = 0.015;
dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;
xxd = 0:XL/(nxd-1):XL;
yyd = 0:YL/(nyd-1):YL;
white = zeros(ny,nx);
br    = zeros(ny,nx); %Radial magnetic flux density [T]
bz    = zeros(ny,nx); %Axial magnetic flux density [T]
babs  = zeros(ny,nx); 
brcen = zeros(1,nx);
brcen2 = zeros(1,nx);
brdata = zeros(nxd*nyd,1);

%The first point is to reproduce Br at centerline

c6 =-1.739302e+8;
c5 = 1.996393e+7;
c4 =-7.588146e+5;
c3 = 8.951234e+3;
c2 = 1.536892e+1;
c1 = 5.324089e-1;
c0 = 3.204751e-3;

for i = 1:1:nx
   if i<=floor(nx*0.70)
      brcen(1,i) = c6*xx(i)^6+c5*xx(i)^5+c4*xx(i)^4+c3*xx(i)^3+c2*xx(i)^2+c1*xx(i)+c0;
   else
      brcen(1,i)=2.5*brcen(1,i-1)-2*brcen(1,i-2)+0.5*brcen(1,i-3);
   end
end

%Shielded magnetic field
bo = nx/2;
for i = bo:1:nx
   brcen2(1,i) = c6*xx(i)^6+c5*xx(i)^5+c4*xx(i)^4+c3*xx(i)^3+c2*xx(i)^2+c1*xx(i)+c0;
end
c6 = 2.27186413E+09;
c5 =-1.33870441E+08;
c4 = 2.37358937E+06;
c3 =-1.18164059E+04;
c2 = 5.51909707E+01;
c1 =-7.66666435E-02;
c0 = 1.99999784E-04;
for i = 1:1:bo-1
   brcen2(1,i) = c6*xx(i)^6+c5*xx(i)^5+c4*xx(i)^4+c3*xx(i)^3+c2*xx(i)^2+c1*xx(i)+c0;
end
%Normally comment out
%brcen(1,:) = brcen2(1,:);

yfunc = zeros(ny,1);
xfunc = zeros(1,nx);
for j = 1:1:ny
   for i = 1:1:nx
      yfunc(j,1) = 1-2*(j-1)/(ny-1);
      xfunc(1,i) = -0.4+1.6*(((i-1)-(nx-1))/(nx-1))^2;
      theta = 0.5*PI+0.4*PI*yfunc(j,1)*xfunc(1,i);
      bstr  = brcen(1,i);
      br(j,i) = bstr*sin(theta);
      bz(j,i) = bstr*cos(theta);
      babs(j,i)= sqrt(bz(j,i)^2+br(j,i)^2);
   end
end


%set(gcf, 'PaperPosition', [2 1 6.8 5.2]);
set(gcf, 'PaperPosition', [2 1 5.1 3.9]);
plot(xx(1,:),brcen(1,:)*1e3,'k-','LineWidth',1.0);
hold on
%plot(xx(1,:),tecen(1,:),'b--','LineWidth',1.0);
%hold on
%plot(xx(1,:),tecen1(1,:),'k--','LineWidth',1.0);
%hold off
set(gca,'XLim', [0,XL]);
set(gca,'YLim', [0,25]);
set(gca,'XTick',[0 0.01 0.02 0.03 0.04 0.05])
set(gca,'XTickLabel',{'0.0','10.0','20.0','30.0','40.0','50.0'})
set(gca,'YTick',[0 5 10 15 20 25])
set(gca,'YTickLabel',{'0.0','5.0','10.0','15.0','20.0','25.0'})
xlabel('Axial position, mm','FontSize',16)
%ylabel({'Magnetic flux density, mT','Electron temperature, eV'},'FontSize',16)
ylabel({'Magnetic flux density, mT'},'FontSize',16)
%legend('Magnetic flux density','Electron temperature','Location','NorthEast')
saveas(figure(1),'/Users/kawashima/Dropbox/HZT2D/figure/datafile/Brcenter.png');
hold off

for j = 1:1:nyd
   for i = 1:1:nxd
      ip = (i-1)*ratio+1;
      jp = (j-1)*ratio+1;
      brdata(nxd*(j-1)+i) = brcen(1,ip);
   end
end
MFdata  = brdata;
dlmwrite(strcat('/Users/kawashima/Dropbox/HZT2D/program/','SPT100MF48.dat'),MFdata, 'delimiter', '\t','precision','%20.10e');


%byy  = reshape(bydata,nxd,nyd)';

