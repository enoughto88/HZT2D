%clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.5 5]);
set(gca,'XTick',[0 0.02 0.04 0.06 0.08 0.1])
set(gca,'XTickLabel',{'0','0.02','0.04','0.06','0.08','0.1'})
%{
ncell =zeros(10,1);
time1 =zeros(10,1);
time2 =zeros(10,1);
ref1  =zeros(10,1);
ref2  =zeros(10,1);
%}
ncell = data(:, 1);
time1 = data(:, 2);
time2 = data(:, 3);
time3 = data(:, 4);
time4 = data(:, 5);

%plot(ncell,time1,'ko',ncell,time2,'k^',ncell,time3,'ks','MarkerSize',7,'LineWidth',1.5);

plot(ncell,time1,'k^',ncell,time2,'ko','MarkerSize',7,'LineWidth',1.5);
hold on
plot(ncell,time3,'k-',ncell(1:5),time4(1:5),'k--','LineWidth',1.0);
%hold on
grid off
set(gca, 'XLim', [2.5, 4.2]);
set(gca, 'YLim', [4.5, 6.5]);
set(gca,'XTick',[2.5 3 3.5 4])
set(gca,'XTickLabel',{'2.5','3.0','3.5','4.0'})
set(gca,'YTick',[4.5 5 5.5 6 6.5])
set(gca,'YTickLabel',{'4.5','5.0','5.5','6.0','6.5'})
xlabel('log_{10}(N_{cell})','FontSize',16)
ylabel('log_{10}(T_{converge})','FontSize',16)
legend('Linearization','LU-ADI','1.5th-order slope','3rd-order slope','Location','NorthWest')
saveas(figure(1),'/Users/kawashima/Desktop/calctime.png');
hold off

