FileToOpen = 'struct_test.dat';

% Read Data
file=fopen(FileToOpen);
C=textscan(file,'%f%f%f%f%f','commentStyle','#');
X = C{1};
Y = C{2};
Z = C{3};
R = C{4};
Phi = C{5};

% Limit Range
X = X(find(Phi>=3/2*pi & Phi < 10*pi));
Y = Y(find(Phi>=3/2*pi & Phi < 10*pi));
Z = Z(find(Phi>=3/2*pi & Phi < 10*pi));

% Plot Wall
wallfile=fopen('C:\c++\d3d\wall.dat');
C3=textscan(wallfile,'%f%f%f%f','commentStyle','#');
XW = C3{3};
YW = zeros(48,1);
ZW = C3{4};

% FontName and FontSize
schrift = 'Times';
schrift_size = 16;

% Plot
figure
axes('FontName',schrift,'FontSize',schrift_size);
plot3(X,Y,Z,'b-','LineWidth',2.5)
hold on
plot3(XW,YW,ZW,'k--','LineWidth',1.5)
plot3(-XW,YW,ZW,'k--','LineWidth',1.5)
plot3(YW,XW,ZW,'k--','LineWidth',1.5)
plot3(YW,-XW,ZW,'k--','LineWidth',1.5)
plot3(X(1),Y(1),Z(1),'rx','LineWidth',1.5,'markerSize',8)
plot3(X(end),Y(end),Z(end),'go','LineWidth',1.5,'markerSize',6)

% Setup Axis
axis tight;
set(gca,'linewidth',1.5,'FontName',schrift);
xlabel({'X [m]'},'FontName',schrift,'FontSize',schrift_size);
ylabel({'Y [m]'},'FontName',schrift,'FontSize',schrift_size);
zlabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
grid on
set(gca,'DataAspectRatio',[1 1 1]);
view(45,20);

% set screen size and position of figure window
hoehe = 13;
set(gcf, 'Units', 'centimeters');
set(0, 'Units', 'centimeters');
scrsz = get(0,'ScreenSize');
set(0, 'Units', 'pixels');
set(gcf, 'Position', [(scrsz(3)-18)/2 (scrsz(4)-15)/2 18 hoehe]);

% Print to file
filenameout = FileToOpen(1:length(FileToOpen)-4);
print('-r300', '-djpeg', filenameout)