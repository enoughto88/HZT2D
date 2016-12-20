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
fileini = 1;
fileend = 2000;
aveini  = 1;
aveend  = 24;
dtpic   = 1.0e-8;
ismp    = 50;
%Grid parameter
nx = 24;
ny = 24;
XL = 50;
YL = 360;
outnion = 1;
maxnion = 2e19;
%----------------------------------


dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;
time = 0:dtpic*ismp:dtpic*ismp*(fileend-fileini);

nneu   = zeros(ny,nx);
nion   = zeros(ny,nx);
Tele   = zeros(ny,nx);
phii   = zeros(ny,nx);
qion   = zeros(ny,nx);
uele   = zeros(ny,nx);
vele   = zeros(ny,nx);
nneuy  = zeros(ny,1);
Teley  = zeros(ny,1);
phiiy  = zeros(ny,1);
qiony  = zeros(ny,1);
nneu_ave   = zeros(ny,nx);
nion_ave   = zeros(ny,nx);
Tele_ave   = zeros(ny,nx);
phii_ave   = zeros(ny,nx);
qion_ave   = zeros(ny,nx);
uele_ave   = zeros(ny,nx);
vele_ave   = zeros(ny,nx);

nion_his   = zeros(nx,fileend-fileini+1);

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

nionx  = zeros(nx,1)';
for jj = aveini:1:aveend
   nionx(1,:) = nionx(1,:)+nion(jj,:)/(aveend-aveini+1);
end
nion_his(:,ii-fileini+1) = nionx(1,:);

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



if outnion == 1
set(gcf, 'PaperPosition', [2 1 13.6 5.2]);
[c,h]=contourf(time*1e3,xx,nion_his/1e18,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Ion number density, 10^{18} m^{-3}','FontSize',16)
shading interp
colormap('jet')
%caxis([0 100])
hold on
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,1.0]);
set(gca, 'YLim', [0,XL]);
set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])
set(gca,'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'})
set(gca,'YTick',[0 10 20 30 40 50])
set(gca,'YTickLabel',{'0','10','20','30','40','50'})
xlabel('Time, {\mu}s','FontSize',16)
ylabel('Axial position, mm','FontSize',16)
saveas(figure(1),strcat(dirout,'nion_his2.png'));
hold off
close all
end







