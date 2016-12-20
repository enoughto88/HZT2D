clear all
close all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5 5.2]);

nx = 24+1;
ny = 24+1;
X0 = 0;
XL = 50;
Y0 = -0.0;
YL = 360.0;
dx = (50-0)/(nx-1);
dy = (360-0)/(ny-1);

dname = '/Users/kawashima/Dropbox/HZT2D';
xx  = zeros(ny,nx);
yy  = zeros(ny,nx);

for j=1:1:ny
    for i=1:1:nx
        xx(j,i) = 0+(i-1)*dx;
        yy(j,i) = 0+(j-1)*dy;
    end
end

for j=1:1:ny
    for i=1:1:nx
        plot(xx(j,:),yy(j,:),'b-')
        hold on
    end
end

for i=1:1:nx
    for j=1:1:ny
        plot(xx(:,i),yy(:,i),'b-')
        hold on
    end
end


xlabel('Axial position, mm','FontSize',16)
ylabel('Azimuthal position, deg','FontSize',16)
set(gca, 'XLim', [X0, XL]);
set(gca, 'YLim', [Y0, YL]);
set(gca,'XTick',[0 10 20 30 40 50])
set(gca,'XTickLabel',{'0','10','20','30','40','50'})
set(gca,'YTick',[0 90 180 270 360])
set(gca,'YTickLabel',{'0','90','180','270','360'})
fname  = strcat(dname,'/figure/datafile/grid.png');
saveas(figure(1),fname);
hold off



