%Particle motion visualization
%Directries are written in Linux format
%Input data file is trajectory_ion.dat or trajectory_ele.dat

%Preambles
clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16)
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.2 5]);

%-----Parameters (Change here)-----
%Directories for input data and output figures
dirin  = '/Users/kawashima/Dropbox/HZT2D/output/datafile/';
dirout = '/Users/kawashima/Dropbox/HZT2D/figure/datafile/';
%Grid parameter
XL = 50;
YL = 360;
%Particle number to display motion (1-50)
pnum = 11;
%Input and output file name
fname  = 'trajectory';
%----------------------------------

%Download particle position data
data=dlmread([strcat(dirin,fname),'.dat']);
leng = length(data(:,1));
x = data(1:leng, 3*pnum-2);
y = data(1:leng, 3*pnum-1);
s = data(1:leng, 3*pnum  );
%Make figure of particle trajectory
plot(x*1e3,y/(2.0d0*pi*0.04d0)*YL,'bo','MarkerSize',2);
grid on
set(gca, 'XLim', [0,XL]);
set(gca, 'YLim', [0,360]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
saveas(figure(1),[strcat(dirout,fname),'.png']);
