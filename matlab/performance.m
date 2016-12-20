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
fname  = 'performance';
%Horizontal axis information
maxnstp  = 2000;
dtpic    = 1.0e-8;
ismp     = 50;
%Calculate average value between avemin and avemax
avemin   = 1;
avemax   = 1;
%Switch for output
outflow  = 1;
outforce = 1;
%Maximum and minimum values for vertical axis
%maxfx = 300.0e-9; minfx =-50.00e-9;
%maxfy = 150.0e-9; minfy =-150.0e-9;
%maxtz = 30.00e-9; mintz =-30.00e-9;
%----------------------------------


%Download performance data
data     = dlmread([strcat(dirin,fname),'.dat']);
nstp     = data(:, 1);
neuflow  = data(:, 2);
ionflow  = data(:, 3);
xforce   = data(:, 4);
acurrent =-data(:, 5);

%Moving average
a = 1;
b = ones(1,10)/10;
neuflowfilt=filter(b,a,neuflow);
neuflowave= zeros(length(data),1);
neuflowave(:,1)= mean(neuflow(avemin:avemax));
ionflowfilt=filter(b,a,ionflow);
ionflowave= zeros(length(data),1);
ionflowave(:,1)= mean(ionflow(avemin:avemax));
xforcefilt=filter(b,a,xforce);
xforceave= zeros(length(data),1);
xforceave(:,1)= mean(xforce(avemin:avemax));
acurrentfilt=filter(b,a,acurrent);
acurrentave= zeros(length(data),1);
acurrentave(:,1)= mean(acurrent(avemin:avemax));


%Make figure of plasma drag force
if(outflow==1)
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.0 3.9]);
plot(nstp*dtpic*ismp*1e3,(ionflowfilt+neuflowfilt)*1e6,'k-')
hold on
plot(nstp*dtpic*ismp*1e3,ionflowfilt*1e6,'b--')
hold on
plot(nstp*dtpic*ismp*1e3,neuflowfilt*1e6,'r-.')
%plot(nstp*dtpic*ismp*1e3,neuflowave*1e9,'k--')
%str1 = ['Average: ' strcat(num2str(floor(neuflowave(1,1)*1e10)/10),' nN')];
%annotation('textbox', [0.58,0.7,0.4,0.2],'String',str1,'FontSize',16,'LineStyle','none')
xlabel('Time, ms','FontSize',16)
ylabel('Mass flow rate, mg s^{-1}','FontSize',16)
set(gca, 'XLim', [0, maxnstp*dtpic*ismp*1e3]);
%set(gca, 'YLim', [minfx, maxfx]*1e9)
legend('Total','Ion','Neutral','Location','NorthEast')
saveas(figure(1),[dirout,'flow.png']);
hold off
end


if(outforce==1)
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.0 3.9]);
plot(nstp*dtpic*ismp*1e3,xforcefilt*1e3,'k-')
%plot(nstp*dtpic*ismp*1e3,neuflowave*1e9,'k--')
%str1 = ['Average: ' strcat(num2str(floor(neuflowave(1,1)*1e10)/10),' nN')];
%annotation('textbox', [0.58,0.7,0.4,0.2],'String',str1,'FontSize',16,'LineStyle','none')
xlabel('Time, ms','FontSize',16)
ylabel('Thrust, mN','FontSize',16)
set(gca, 'XLim', [0, maxnstp*dtpic*ismp*1e3]);
set(gca, 'YLim', [0, 150])
%legend('Total','Ion','Neutral','Location','NorthEast')
saveas(figure(1),[dirout,'thrust.png']);
hold off
end
%{
if(outtz==1)
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.5 4.5]);
saveas(figure(1),[dirout,'torquez.png']);
plot(nstp*dtpic*ismp*1e3,ztorquefilt*1e9,'b-')
hold on
plot(nstp*dtpic*ismp*1e3,ztorqueave*1e9,'k--')
str2 = ['Average: ' strcat(num2str(floor(ztorqueave(1,1)*1e10)/10),' nNm')];
annotation('textbox', [0.54,0.7,0.4,0.2],'String',str2,'FontSize',16,'LineStyle','none')
xlabel('Time, ms','FontSize',16)
ylabel('z-Torque, nNm','FontSize',16)
set(gca, 'XLim', [0, maxnstp*dtpic*ismp*1e3]);
set(gca, 'YLim', [mintz, maxtz]*1e9)
%set(gca,'XTick',[0.0 2.0 4.0 6.0 8.0 10.0 12.0])
%set(gca,'XTickLabel',{'0.0','2.0','4.0','6.0','8.0','10.0','12.0'})
%set(gca,'YTick',[0 5 10 15 20 25])
%set(gca,'YTickLabel',{'0.0','5.0','10.0','15.0','20.0','25.0'})
saveas(figure(1),[dirout,'torquez.png']);
hold off
close all
end
%}





