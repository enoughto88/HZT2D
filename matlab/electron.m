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
fname1 = 'electron';
%Average the data of fileini ~ fileend
fileini = 1;
fileend = 2000;
%Grid parameter
nx = 24;
ny = 24;
XL = 50;
YL = 360;
outuele = 0;
outcond = 1;
outOmee = 1;
outfion = 1;
outfcol = 1;
outdist = 0;
%----------------------------------


dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;

uele   = zeros(ny,nx);
cond   = zeros(ny,nx);
Omee   = zeros(ny,nx);
fion   = zeros(ny,nx);
fcol   = zeros(ny,nx);
uelec  = zeros(1,nx);
condc  = zeros(1,nx);
Omeec  = zeros(1,nx);
fionc  = zeros(1,nx);
fcolc  = zeros(1,nx);
uele_ave   = zeros(ny,nx);
cond_ave   = zeros(ny,nx);
Omee_ave   = zeros(ny,nx);
fion_ave   = zeros(ny,nx);
fcol_ave   = zeros(ny,nx);

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

uele  = reshape(data(:,1),nx,ny)';
cond  = reshape(data(:,2),nx,ny)';
Omee  = reshape(data(:,3),nx,ny)';
fion  = reshape(data(:,4),nx,ny)';
fcol  = reshape(data(:,5),nx,ny)';

uele_ave = uele_ave+uele;
cond_ave = cond_ave+cond;
Omee_ave = Omee_ave+Omee;
fion_ave = fion_ave+fion;
fcol_ave = fcol_ave+fcol;
end
uele = uele_ave/(fileend-fileini+1);
cond = cond_ave/(fileend-fileini+1);
Omee = Omee_ave/(fileend-fileini+1);
fion = fion_ave/(fileend-fileini+1);
fcol = fcol_ave/(fileend-fileini+1);


if outcond == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,cond/1e23,20);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Electron conductivity, 10^{23}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
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
saveas(figure(1),strcat(dirout,'cond.png'));
hold off
close all
end

if outOmee == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,Omee/1e4,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Electron Hall parameter, 10^4','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
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
saveas(figure(1),strcat(dirout,'Omee.png'));
hold off
close all
end

if outfion == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,fion/1e5,20);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Ionization collision frequency, 10^{5} s^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
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
saveas(figure(1),strcat(dirout,'fion.png'));
hold off
close all
end
if outfcol == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,fcol/1e7,20);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Total collision frequency, 10^{7} s^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 300])
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
saveas(figure(1),strcat(dirout,'fcol.png'));
hold off
close all
end

%{
if outdist == 1
for j = 1:1:ny
for i = 1:1:nx
   nneuc(i)  = nneuc(i)+uele(j,i)/ny;
   nionc(i)  = nionc(i)+vele(j,i)/ny;
   Telec(i)  = Telec(i)+Omee(j,i)/ny;
   phiic(i)  = phiic(i)+fion(j,i)/ny;
   qionc(i)  = qionc(i)+fcol(j,i)/ny;
end
end
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.5 4.5]);
%plot(xx,nneuc/max(nneuc),'k-')
%hold on
plot(xx,Telec/max(Telec),'k-')
hold on
plot(xx,nionc/max(nionc),'b--')
%hold on
%plot(xx,phiic/max(phiic),'k-')
hold on
plot(xx,qionc/max(qionc),'r-.')
hold on
xlabel('Axial position, mm','FontSize',16)
ylabel('Norm. plasma property','FontSize',16)
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,1.1]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})
legend('Te','n_e','q_{ion}','NorthEast')
plot([23 23], [0 1.1], 'k--')
saveas(figure(1),strcat(dirout,'dist.png'));
hold off
close all
end
%}

