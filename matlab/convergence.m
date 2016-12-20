clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5 4]);
set(gca,'XTick',[0 0.02 0.04 0.06 0.08 0.1])
set(gca,'XTickLabel',{'0','0.02','0.04','0.06','0.08','0.1'})

fileini = 310;
nfile   = 310;

fname = ['/Users/kawashima/Dropbox/HZT2D/output/datafile/convergence'];

data = dlmread([fname,'.dat']);
nstp = data(:, 1);
res1 = data(:, 2);
res2 = data(:, 3);
res3 = data(:, 4);


semilogy(nstp,res1,'k-',nstp,res2,'b--',nstp,res3,'r-.','LineWidth',1.5);
grid off
set(gca, 'XLim', [0, 150000]);
set(gca, 'YLim', [1e-20, 1e0]);
xlabel('Number of time steps','FontSize',16)
ylabel('Normalized difference, D_{norm}','FontSize',16)
legend('Space potential','z-Momentum','\theta-Momentum')
saveas(figure(1),'/Users/kawashima/Dropbox/HZT2D/figure/datafile/convergence.png');