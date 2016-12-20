%Plasma force and torque visualization
%Directries are written in Linux format
%Input data file is performance.dat

%Preambles
clear all
close all

%-----Parameters (Change here)-----
%Directories for input data and output figures
dirin  = '/Users/kawashima/Dropbox/HZT2D/output/datafile/';
dirout = '/Users/kawashima/Dropbox/HZT2D/figure/datafile/';
fname  = 'resn';
%Horizontal axis information
maxnstp  = 1000;
dtpic    = 1.0e-8;
ismp     = 50;
%Calculate average value between avemin and avemax
avemin   = 1;
avemax   = 1;
%Switch for output
outresn = 1;
%Maximum and minimum values for vertical axis
%maxfx = 300.0e-9; minfx =-50.00e-9;
%maxfy = 150.0e-9; minfy =-150.0e-9;
%maxtz = 30.00e-9; mintz =-30.00e-9;
%----------------------------------


%Download performance data
data     = dlmread([strcat(dirin,fname),'.dat']);
nstp     = data(:, 1);
lpele  = data(:, 2);
resn1  = data(:, 3);
resn2  = data(:, 4);
resn3  = data(:, 5);

%Moving average
a = 1;
b = ones(1,10)/10;
resn1filt=filter(b,a,resn1);
resn2filt=filter(b,a,resn2);
resn3filt=filter(b,a,resn3);
%resn1ave= zeros(length(data),1);
%resn1ave(:,1)= mean(resn1(avemin:avemax));

%Make figure of plasma drag force
if(outresn==1)
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 6.8 5.2]);
semilogy(nstp*dtpic*ismp*1e3,resn1filt,'k-')
hold on
semilogy(nstp*dtpic*ismp*1e3,resn2filt,'b--')
hold on
semilogy(nstp*dtpic*ismp*1e3,resn3filt,'r-.')
%plot(nstp*dtpic*ismp*1e3,neuflowave*1e9,'k--')
%str1 = ['Average: ' strcat(num2str(floor(neuflowave(1,1)*1e10)/10),' nN')];
%annotation('textbox', [0.58,0.7,0.4,0.2],'String',str1,'FontSize',16,'LineStyle','none')
xlabel('Time, ms','FontSize',16)
ylabel('Normalized difference of variable  {\it D}_{norm}','FontSize',16)
set(gca, 'XLim', [0, maxnstp*dtpic*ismp*1e3]);
set(gca, 'YLim', [1e-8, 1e-4])
legend('Space potential','x-Electron momentum','y-Electron momentum','Location','NorthEast')
saveas(figure(1),[dirout,'resn.png']);
hold off
end







