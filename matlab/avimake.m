close all
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
dirout = '/Users/kawashima/Desktop/movie/';
fname1 = 'distribution';
%Average the data of fileini ~ fileend
fileini = 561;
fileend = 760;
dtpic   = 1.0e-8;
ismp    = 50;
%Grid parameter
nx = 24;
ny = 24;
XL = 50;
YL = 360;
maxnion = 15e18;
maxqion = 1.5e24;
maxphi  = 310;
maxnneu = 8e19;
pic = 0;
movie = 1;
nionpic = 0;
qionpic = 0;
phipic  = 0;
nneupic = 0;
%----------------------------------


dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;
xx2= 0:XL/(nx-3):XL;


nion1  = zeros(ny,nx);
nion2  = zeros(ny,nx);
nion3  = zeros(ny,nx);
nion4  = zeros(ny,nx);
nion   = zeros(ny,nx);
qion1  = zeros(ny,nx);
qion2  = zeros(ny,nx);
qion3  = zeros(ny,nx);
qion4  = zeros(ny,nx);
qion   = zeros(ny,nx);
phi1  = zeros(ny,nx);
phi2  = zeros(ny,nx);
phi3  = zeros(ny,nx);
phi4  = zeros(ny,nx);
phi   = zeros(ny,nx);
nneu1  = zeros(ny,nx);
nneu2  = zeros(ny,nx);
nneu3  = zeros(ny,nx);
nneu4  = zeros(ny,nx);
nneu   = zeros(ny,nx);

if pic==1
for ii=fileini:2:fileend
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
nion1  = reshape(data(:,2),nx,ny)';
qion1  = reshape(data(:,5),nx,ny)';
phi1   = reshape(data(:,4),nx,ny)';
nneu1  = reshape(data(:,1),nx,ny)';

if ii+1<10
   fname2  = '.00000';
elseif ii+1<100
   fname2  = '.0000';
elseif ii+1<1000
   fname2  = '.000';
elseif ii+1<10000
   fname2  = '.00';
else
   fname2  = '.0';
end
fname = strcat(dirin,fname1,fname2,num2str(ii+1));
data  = dlmread([fname,'.dat']);
nion2  = reshape(data(:,2),nx,ny)';
qion2  = reshape(data(:,5),nx,ny)';
phi2   = reshape(data(:,4),nx,ny)';
nneu2  = reshape(data(:,1),nx,ny)';

if ii+2<10
   fname2  = '.00000';
elseif ii+2<100
   fname2  = '.0000';
elseif ii+2<1000
   fname2  = '.000';
elseif ii+2<10000
   fname2  = '.00';
else
   fname2  = '.0';
end
fname = strcat(dirin,fname1,fname2,num2str(ii+2));
data  = dlmread([fname,'.dat']);
nion3  = reshape(data(:,2),nx,ny)';
qion3  = reshape(data(:,5),nx,ny)';
phi3   = reshape(data(:,4),nx,ny)';
nneu3  = reshape(data(:,1),nx,ny)';

if ii+3<10
   fname2  = '.00000';
elseif ii+3<100
   fname2  = '.0000';
elseif ii+3<1000
   fname2  = '.000';
elseif ii+3<10000
   fname2  = '.00';
else
   fname2  = '.0';
end
fname = strcat(dirin,fname1,fname2,num2str(ii+3));
data  = dlmread([fname,'.dat']);
nion4  = reshape(data(:,2),nx,ny)';
qion4  = reshape(data(:,5),nx,ny)';
phi4   = reshape(data(:,4),nx,ny)';
nneu4  = reshape(data(:,1),nx,ny)';

nion = 0.25*(nion1+nion2+nion3+nion4);
qion = 0.25*(qion1+qion2+qion3+qion4);
phi  = 0.25*(phi1+phi2+phi3+phi4);
nneu = 0.25*(nneu1+nneu2+nneu3+nneu4);

if nionpic == 1
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
caxis([0 maxnion/1e18])
set(h,'YTick',[0 5 10 15 20 25 30])
set(h,'YTickLabel',{'0','5','10','15','20','25','30'})
hold on
str2 = ['t = {       } \mus'];
str3 = [strcat(num2str(floor(ii*dtpic*ismp*1e6+0.5)+1))];
if ii==1
   str3 = [strcat(num2str(floor(0*dtpic*ismp*1e6)))];
end
if ii==1995
   str3 = [strcat(num2str(floor(2000*dtpic*ismp*1e6)))];
end
annotation('textbox', [0.58,0.72,0.4,0.2],'String',str2,'FontSize',16,'Color','w','LineStyle','none')
hold on
annotation('textbox', [0.63,0.72,0.4,0.2],'String',str3,'FontSize',16,'Color','w','LineStyle','none')
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
fname = strcat(dirout,'mov_',num2str(ii));
saveas(figure(1),strcat(fname,'.png'));
hold off
close all
end

if qionpic == 1
qion(:,:) = min(qion(:,:),maxqion);
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,qion/1e23,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'Ion production rate, 10^{23} m^{-3} s^{-1}','FontSize',16)
shading interp
colormap('jet')
caxis([0 maxqion/1e23])
set(h,'YTick',[0 5 10 15])
set(h,'YTickLabel',{'0','5','10','15'})
hold on
str2 = ['t = {       } \mus'];
str3 = [strcat(num2str(floor(ii*dtpic*ismp*1e6+0.5)+1))];
if ii==1995
   str3 = [strcat(num2str(floor(2000*dtpic*ismp*1e6)))];
end
annotation('textbox', [0.58,0.72,0.4,0.2],'String',str2,'FontSize',16,'Color','w','LineStyle','none')
hold on
annotation('textbox', [0.63,0.72,0.4,0.2],'String',str3,'FontSize',16,'Color','w','LineStyle','none')
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
fname = strcat(dirout,'mov_',num2str(ii));
saveas(figure(1),strcat(fname,'.png'));
hold off
close all
end

if phipic == 1
phi(:,:) = min(phi(:,:),maxphi);
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,phi,100);
set(h,'edgecolor','none')
grid off
h = colorbar;
ylabel(h,'Space potential, V','FontSize',16)
shading interp
colormap('jet')
caxis([0 maxphi])
set(h,'YTick',[0 100 200 300])
set(h,'YTickLabel',{'0','100','200','300'})
hold on
str2 = ['t = {       } \mus'];
str3 = [strcat(num2str(floor(ii*dtpic*ismp*1e6+0.7)+1))];
if ii==1
   str3 = [strcat(num2str(floor(0*dtpic*ismp*1e6)))];
end
if ii==1995
   str3 = [strcat(num2str(floor(2000*dtpic*ismp*1e6)))];
end
annotation('textbox', [0.58,0.72,0.4,0.2],'String',str2,'FontSize',16,'Color','w','LineStyle','none')
hold on
annotation('textbox', [0.63,0.72,0.4,0.2],'String',str3,'FontSize',16,'Color','w','LineStyle','none')
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
fname = strcat(dirout,'mov_',num2str(ii));
saveas(figure(1),strcat(fname,'.png'));
hold off
close all
end
if nneupic == 1
nneu(:,:) = min(nneu(:,:),maxnneu);
set(gcf, 'PaperPosition', [2 1 5.85 5.2]);
[c,h]=contourf(xx,yy,nneu/1e19,15);
set(h,'edgecolor','k')
grid off
h = colorbar;
%ylabel(h,'\phi / \phi_{Anode}','FontSize',16)
ylabel(h,'Neutral number density, 10^{19} m^{-3}','FontSize',16)
%shading interp
colormap('jet')
caxis([0 maxnneu/1e19])
set(h,'YTick',[0 5 10 15 20 25 30])
set(h,'YTickLabel',{'0','5','10','15','20','25','30'})
hold on
str2 = ['t = {       } \mus'];
str3 = [strcat(num2str(floor(ii*dtpic*ismp*1e6+0.5)+1))];
if ii==1
   str3 = [strcat(num2str(floor(0*dtpic*ismp*1e6)))];
end
if ii==1995
   str3 = [strcat(num2str(floor(2000*dtpic*ismp*1e6)))];
end
annotation('textbox', [0.58,0.72,0.4,0.2],'String',str2,'FontSize',16,'Color','w','LineStyle','none')
hold on
annotation('textbox', [0.63,0.72,0.4,0.2],'String',str3,'FontSize',16,'Color','w','LineStyle','none')
%plot([23 23], [0 YL], 'w--')
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
fname = strcat(dirout,'mov_',num2str(ii));
saveas(figure(1),strcat(fname,'.png'));
hold off
close all
end
end
end

if movie == 1
for ii=fileini:2:fileend
   fname = strcat(dirout,'mov_',num2str(ii));
   fname = strcat(fname,'.png');
   moviepic=imread(fname);
   fnum = (ii-fileini)/2+1;
   mov(fnum) =im2frame((moviepic));
end
movie2avi(mov,strcat(dirout,'movie.avi'),'fps',25,'compression','None','quality',100)
end





