clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');


%-----Parameters (Change here)-----
%Directories for input data and output figures
dirin  = '/Users/kawashima/Dropbox/HZT2D/output/distribution/';
dirout = '/Users/kawashima/Dropbox/HZT2D/figure/distribution/';
fname1 = 'distribution';
%Average the data of fileini ~ fileend
fileini = 1804;
fileend = 1805;
dtpic   = 1.0e-8;
ismp    = 50;
%Grid parameter
nx = 24;
ny = 24;
XL = 50;
YL = 360;
outnneu = 1;
outnion = 1;
outTele = 1;
outphii = 1;
outqion = 1;
outdist = 0;
outuele = 1;
outflow = 1;
maxnion = 1e20;
%----------------------------------


dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;


nneu   = zeros(ny,nx);
nion   = zeros(ny,nx);
Tele   = zeros(ny,nx);
phii   = zeros(ny,nx);
qion   = zeros(ny,nx);
uele   = zeros(ny,nx);
vele   = zeros(ny,nx);
nneuc  = zeros(1,nx);
nionc  = zeros(1,nx);
Telec  = zeros(1,nx);
phiic  = zeros(1,nx);
qionc  = zeros(1,nx);
nneu_ave   = zeros(ny,nx);
nion_ave   = zeros(ny,nx);
Tele_ave   = zeros(ny,nx);
phii_ave   = zeros(ny,nx);
qion_ave   = zeros(ny,nx);
uele_ave   = zeros(ny,nx);
vele_ave   = zeros(ny,nx);

for ii=fileini:1:fileend
if ii<10
   fname2  = '.00000';
elseif ii<100
   fname2  = '.0000';
elseif ii<1000
   fname2  = '.000';
elseif ii<10000
   fname2  = '.00';
else
   fname2  = '.0';
end

fname = strcat(dirin,fname1,fname2,num2str(ii));
data  = dlmread([fname,'.dat']);

nneu  = reshape(data(:,1),nx,ny)';
nion  = reshape(data(:,2),nx,ny)';
Tele  = reshape(data(:,3),nx,ny)';
phii  = reshape(data(:,4),nx,ny)';
qion  = reshape(data(:,5),nx,ny)';
uele  = reshape(data(:,6),nx,ny)';
vele  = reshape(data(:,7),nx,ny)';

nneu_ave = nneu_ave+nneu;
nion_ave = nion_ave+nion;
Tele_ave = Tele_ave+Tele;
phii_ave = phii_ave+phii;
qion_ave = qion_ave+qion;
uele_ave = uele_ave+uele;
vele_ave = vele_ave+vele;
end
nneu = nneu_ave/(fileend-fileini+1);
nion = nion_ave/(fileend-fileini+1);
Tele = Tele_ave/(fileend-fileini+1);
phii = phii_ave/(fileend-fileini+1);
qion = qion_ave/(fileend-fileini+1);
uele = uele_ave/(fileend-fileini+1);
vele = vele_ave/(fileend-fileini+1);


if outnneu == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,nneu/1e19,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'Neutral number density, 10^{19} m^{-3}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 100])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'nneu.png'));
hold off
close all
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,log10(nneu),20);
set(h,'edgecolor','none')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'log Neutral number density, 10^{19} m^{-3}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 100])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'lognneu.png'));
hold off
close all
end

if outnion == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
logne = log10(nion);
logne(:,:) = min(logne(:,:),19);
[c,h]=contourf(xx,yy,logne,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'log_{10} (n_e)','FontSize',16)
shading flat
colormap('jet')
caxis([16 max(max(logne))])
set(h,'YTick',[16 17 18 19 20])
set(h,'YTickLabel',{'16','17','18','19','20'})
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'logne.png'));
hold off
nion(:,:) = min(nion(:,:),maxnion);
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,nion/1e18,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'Ion number density, 10^{18} m^{-3}','FontSize',16)
shading interp
colormap('jet')
caxis([0 17])
set(h,'YTick',[0 5 10 15 20])
set(h,'YTickLabel',{'0','5','10','15','20'})
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'nion.png'));
hold off
close all
end

if outTele == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,Tele,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Electron temperature, eV','FontSize',16)
%shading interp
colormap('jet')
caxis([0 max(max(Tele))])
hold on
plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'Tele.png'));
hold off
close all
end


if outphii == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,phii,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Space potential, V','FontSize',16)
shading interp
colormap('jet')
caxis([0 max(max(phii))])
%hold on
%plot([23 23], [0 YL], 'w--')
set(h,'YTick',[0 100 200 300])
set(h,'YTickLabel',{'0','100','200','300'})
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'phii.png'));
hold off
close all
end



if outqion == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,log10(qion),100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'log_{10} (q_{ion})','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'logqion.png'));
hold off

set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,qion/1e23,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Ion production rate, 10^{23} m^{-3} s^{-1}','FontSize',16)
%shading interp
colormap('jet')
caxis([0 15])
set(h,'YTick',[0 2 4 6 8 10 12])
set(h,'YTickLabel',{'0','2','4','6','8','10','12'})
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'qion.png'));
hold off
close all
end


if outuele == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,uele/1e8,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'x-Electron velocity, 10^{8} m s^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'uele.png'));
hold off
close all
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,vele/1e8,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'y-Electron velocity, 10^{8} m s^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'vele.png'));
hold off
close all
end

if outflow == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,phii,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Space potential, V','FontSize',16)
shading interp
colormap('jet')
caxis([0 max(max(phii))])
%hold on
%plot([23 23], [0 YL], 'w--')
set(h,'YTick',[0 100 200 300])
set(h,'YTickLabel',{'0','100','200','300'})
hold on
nline = 10;
ystart=zeros(1,nline)+YL/2;
xstart=[XL/nline/2:XL/nline:XL-XL/nline/2];
h = streamline(xx,yy,uele,vele,xstart,ystart);
set(h,'LineStyle','-','Color','w','LineWidth',0.3)
hold on
h = streamline(xx,yy,-uele,-vele,xstart,ystart);
set(h,'LineStyle','-','Color','w','LineWidth',0.3)
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'flow.png'));
hold off
close all
end


if outdist == 1
for j = 1:1:ny
for i = 1:1:nx
   nneuc(i)  = nneuc(i)+nneu(j,i)/ny;
   nionc(i)  = nionc(i)+nion(j,i)/ny;
   Telec(i)  = Telec(i)+Tele(j,i)/ny;
   phiic(i)  = phiic(i)+phii(j,i)/ny;
   qionc(i)  = qionc(i)+qion(j,i)/ny;
end
end
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 6.8 5.2]);
plot(xx,phiic/max(phiic),'k-')
hold on
%plot(xx,Telec/max(Telec),'k-')
%hold on
plot(xx,nionc/max(nionc),'b--')
hold on
%plot(xx,nneuc/max(nneuc),'g-.')
%hold on
plot(xx,qionc/max(qionc),'r-.')
%hold on
xlabel('Axial position, mm','FontSize',16)
ylabel('Norm. plasma property','FontSize',16)
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,1.1]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})
legend('Space potential','Ion number density','Ion production rate','NorthEast')
%plot([23 23], [0 1.1], 'k--')
saveas(figure(1),strcat(dirout,'dist.png'));
hold off
close all
end


