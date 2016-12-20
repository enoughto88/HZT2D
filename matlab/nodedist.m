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
fname1 = 'nodedist';
%Input file number
file = 5;
%Grid parameter
nx = 24+1;
ny = 24+1;
XL = 50;
YL = 360;
outefnd = 1;
%----------------------------------


dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;

efx   = zeros(ny,nx);
efy   = zeros(ny,nx);

ii = file;
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

efx  = reshape(data(:,1),nx,ny)';
efy  = reshape(data(:,2),nx,ny)';


if outefnd == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,efx,20);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'x-Electric field, V m^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 100])
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'efx.png'));
hold off
close all
end

if outefnd == 1
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,efy,20);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'y-Electric field, V m^{-1}','FontSize',16)
%shading interp
colormap('jet')
%caxis([0 100])
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),strcat(dirout,'efy.png'));
hold off
close all
end
