clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');

dirout = '/Users/kawashima/Dropbox/HZT2D/figure/distribution/';
pi = 3.1415926535893;
nx = 240;
ny = 240;
XL = 50;
YL = 360;
%dx = XL/nx;
%dy = YL/ny;
%xx = dx/2:dx:XL-dx/2;
%yy = dy/2:dy:YL-dy/2;
dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;
ne  = zeros(ny,nx);
te  = zeros(ny,nx);

for j = 1:1:ny
   for i = 1:1:nx
      ne(j,i) = 7e18*exp(-(8*(xx(i)-0.2*XL)/XL)^2)*(1.1*cos(2*pi*((yy(j)-0.5*YL)/YL))+1.2);
      %nn(j,i) = 13e19*(exp(-(5*(xx(i)-0.0*XL)/XL)^2)+0.1)-3e18*exp(-(8*(xx(i)-0.3*XL)/XL)^2)*(1.1*cos(2*pi*((yy(j)-0.5*YL)/YL))+1.2);
      shift =  0.0;
      if xx(i)/XL<shift
        nn(j,i) = 15e19*1.1;
      else
        nn(j,i) = 15e19*(exp(-(3.8*(xx(i)-shift*XL)/XL)^2)+0.1);
      end
   end
end

pos = 1.0;
Tep = 20.0;
Tem = 4.0;
sig = 0.20;
for j = 1:1:ny
   for i = 1:1:nx
      te(j,i) = (Tep-Tem)*exp(-1*(xx(i)*2/XL-pos)^2/(2*sig^2))+Tem;
   end
end
pos = 1.0;
Tep = 20.0;
Tem = 4.0;
sig = 0.20;
for j = 1:1:ny
   for i = 1:1:nx
      te(j,i) = (Tep-Tem)*exp(-1*(xx(i)*2/XL-pos)^2/(2*sig^2))+Tem;
   end
end

set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,ne/1e18,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Ion number density, 10^{18} m^{-3}','FontSize',16)
shading interp
colormap('jet')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'inputne.png'));
hold off


set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,nn/1e19,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Neutral number density, 10^{19} m^{-3}','FontSize',16)
shading interp
colormap('jet')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'inputnn.png'));
hold off

set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,te,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Electron temperature, eV','FontSize',16)
shading interp
colormap('jet')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'inputTe.png'));
hold off

set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,te,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Electron temperature, eV','FontSize',16)
shading interp
colormap('jet')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'inputfion.png'));
hold off

