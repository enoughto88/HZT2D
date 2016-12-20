%Initial condition of phi and te
%
clear all
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [5 5 8.3 3.7]);


PI = 3.1415926535893;
nx = 240;
ny = 240;
nxd= 24;
nyd= 24;
ratio = 10;
XL = 0.05;
YL = 0.015;
dx = XL/(nx-1);
dy = YL/(ny-1);
xx = 0:dx:XL;
yy = 0:dy:YL;
phi   = zeros(ny,nx); %Initial space potential [V]
te    = zeros(ny,nx); %Initial electron temperature [eV]
vy    = zeros(ny,nx); %Radial magnetic flux density [T]
vx    = zeros(ny,nx); %Axial magnetic flux density [T]
brcen = zeros(1,nx);
phicen= zeros(1,nx);
tecen = zeros(1,nx);
pos   = zeros(ny,nx);
phdata = zeros(nxd*nyd,1);
tedata = zeros(nxd*nyd,1);

%The first point is to reproduce Br at centerline
c6 =-1.739302e+8;
c5 = 1.996393e+7;
c4 =-7.588146e+5;
c3 = 8.951234e+3;
c2 = 1.536892e+1;
c1 = 5.324089e-1;
c0 = 3.204751e-3;
for i = 1:1:nx
   if i<=floor(nx*0.70)
      brcen(1,i) = c6*xx(i)^6+c5*xx(i)^5+c4*xx(i)^4+c3*xx(i)^3+c2*xx(i)^2+c1*xx(i)+c0;
   else
      brcen(1,i)=2.5*brcen(1,i-1)-2*brcen(1,i-2)+0.5*brcen(1,i-3);
   end
end

xfuncp = zeros(1,nx);
for i = 1:1:nx
   xfuncp(1,i)  = ((i-1)/(nx-1)).^3;
end
intsum = sum(brcen(1,1:nx-1).^6.*xfuncp(1,1:nx-1));
phicen(1,1) = 300;
for i = 2:1:nx
   phicen(1,i) = phicen(1,i-1)-300*brcen(1,i-1).^6.*xfuncp(1,i-1)/intsum;
end


xfunct = zeros(1,nx);
for i = 1:1:nx
   xfunct(1,i)  = ((i-1)/(nx-1)).^0;
end
intsum = sum(brcen(1,i).^3*xfunct(1,:));
shift = floor(nx/20)+0;
for i = shift:1:nx
   tecen(1,i) =28*brcen(1,i-shift+1).^3*xfunct(1,i-shift+1)/intsum+2;
end
for i = 1:1:shift-1
   tecen(1,i) = tecen(1,shift);
end


for j = 1:1:nyd
   for i = 1:1:nxd
      ip = (i-1)*ratio+1;
      jp = (j-1)*ratio+1;
      phdata(nxd*(j-1)+i) = phicen(1,ip);
      tedata(nxd*(j-1)+i) = tecen(1,ip);
   end
end
MFdata  = horzcat(phdata,tedata);
dlmwrite(strcat('/Users/kawashima/Dropbox/HZT2D/program/','initialphite.dat'),MFdata, 'delimiter', '\t','precision','%20.10e');





