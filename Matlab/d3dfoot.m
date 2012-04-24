clear;
%---------------- Set parameter here -------------------------------
% Color scaling in b
TypeOfPlot = 1;   % 0: contourf  1: pcolor (fast)
autoscale = 2;    %0: use b    1: scale color automatically   2: use b and set SOL value in Colormap to white
b = 0.075:0.005:0.4;   % (con. length only) [0.075:0.005:0.4]  [0.028:0.005:0.4]  [0.028:0.0001:0.108]

% Main Filename
Laminar = 1;          % 0: Footprint  1: Laminar Plot
WhatShallIPlot = 1;   % 0: ntor  1: con. length[km]  2: psimin
printme = 0;          % 0: no export to jpg file     1: export to jpg
FileToOpen = 'lam_test.dat';

% (printme==1 only) Comment line if not wanted: Add an additional string to output filename
%praefix = 'a';

% (Laminar == 0 only) Specify Target and Axis handling 
Target = 1;   % 0: CP  1: Inner + CP  2: Outer  3: Shelf
Physical = 1; % 1: y-Axis in R or Z, depending on target    0: t-Axis   2: t-Axis in cm (inner target only)   
Machine = 1;  % 0: x-Axis in rad 1: x-Axis in deg and REVERSED! (machine coordinates)
CameraView = 0; % 0: no Modification 1: Aspect ratio is set, y-axis covers all tiles

% (Laminar == 0 only) Phi angle area 
VaryPhiRange = 0;  % 0: no Variation 1: Variation 
phimin = (360-105)/180*pi;  %-91
phimax = (360-95)/180*pi; %-43

% Plot additional points 0: none, >=1: Number of Files to plot
plot_manifold = 0;
filenames = {'current_tube5.dat' ...   out   tiles_at_60deg_all.dat   man_unstr1_C_2kA_0_RZ.dat
    'current_tube2.dat' ...         in
    'current_tube3.dat' ...        out
    'current_tube4.dat' ...        out
    'current_tube6.dat' ...        out
    'current_7.dat' ...        out
    'current_tube8.dat'}; ...       in

% Set style for additional data (see Help: LineSpec)
styles = {'w.', 'ko', 'mx', 'm.', 'wx', 'w.', 'k+'};
groesse = [15, 6, 10, 15, 10, 15, 6];    % markersizes for additional data 

% (Laminar==1 & printme==1 only) Adjust horizontal colorbar position in [cm], if necessary
correct_pos = 0.0;    % Default: = 0

    
%------------------------------------------------------------------
%------------------------------------------------------------------
% Code draws color-contour plot of footprint and laminar data
% !!! IMPORTANT !!!
% In footprints: x-Axis (toroidal angle) has to be reversed!
% Code calculates in a right-handed coordiante system 
% -> tor. angle increases counter-clock-wise
% In DIII-D: tor. angle is left-handed -> increases clock-wise
% e.g. 60° in the machine are -60° in footprint and vice versa
%
% if(printme == 1)
%     Figure appearance is prepared for export to file
%     Automatic export to jpg
%     
%     Colorbar position in Laminar plot not perfectly aligned 
%     (reason unknown), but optimized for typical aspect ratio
%     Manual correction afterwards:
%       set(cbar,'location','manual','Position', [11.7+gcapos(1)+<value>, cbpos(2), 0.5, cbpos(4)]);
%       print('-r300', '-djpeg', filenameout)
%     
% end
%
% Warning: in laminar plots check if plot_manifold is either =0
% or only valid data-files are used, otherwise fatal error possible
% (reason unknown!)
%------------------------------------------------------------------
%------------------------------------------------------------------

% Specify Data Format
yVariedFirst = 0;     % Default: =0     use =1 only if file format varies y values first

% Specify b for other than con. length
if(WhatShallIPlot == 2), b = 0.9:0.001:1;  % [0.9:0.001:1]
elseif(WhatShallIPlot == 0) 
    b = 4:0.5:40;
    if(autoscale == 2), autoscale = 0; end
end

% plot limit line between CP and inner target: 0=no 1=yes
if(Laminar==0 && Target==1), plot_target_limit = 1;
else plot_target_limit = 0;
end

% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Check for possible file replacement
if(printme==1)
    if(exist('praefix','var')==0), praefix = ''; end
    filenameout = [FileToOpen(1:length(FileToOpen)-4) praefix];
    if(exist([filenameout '.jpg'],'file')==2)
        disp('File already exists. Enter additional string for Filename.')
        disp('Press just Enter to replace File');
        reply = input('String: ','s');
        if isempty(reply)
            reply = '';
        else filenameout = [filenameout '_' reply];
        end
    end
end

% --- Read data from input file -----------------------------------
file=fopen(FileToOpen);
% die Matrizen hier sind transponiert zu den unteren, ist aber egal
%if(WhatShallIPlot==2)
    C=textscan(file,'%f%f%f%f%f','commentStyle','#');
%else
%    C=textscan(file,'%f%f%f%f%*[^\n]','commentStyle','#');
%end
fclose(file);
phi=C{1};
t=C{2};
val=C{WhatShallIPlot+3};
sizeofvec=length(t);

% transform t to physical values: length in m
if(Laminar==0 && Physical>=1)
    if(Target==1 && Physical == 1) 
        tg = t(t>=0);
        tl = t(t<0);
        tg = -1.22884 - tg*0.13756;
        tl = -1.22884 - tl*0.193967;
        t = [tl;tg];
    elseif(Target==1 && Physical == 2), t = t*19.3967;
    elseif(Target==0), t = -1.22884 - t*0.193967;
    elseif(Target==2), t = 1.15285 + t*0.21915;
    elseif(Target==3), t = 1.372 + t*0.21915;
    end
end

% # grid points of phi
Np=1;
for i=2:sizeofvec
    if(t(i)==t(i-1)) 
        Np=Np+1;
    else break;
    end
end
Nt=sizeofvec/Np;    % # grid points of t

if(yVariedFirst==1)
    Nt=1;
    for i=2:sizeofvec
        if(phi(i)==phi(i-1)) 
            Nt=Nt+1;
        else break;
        end
    end
    Np=sizeofvec/Nt;    % reset # grid points of phi
end

% Set 2D color data
X = reshape(phi,Np,Nt);
Y = reshape(t,Np,Nt);
Z = reshape(val,Np,Nt);

% --- change Phi range --------------------------------------------
if(Laminar==0 && VaryPhiRange==1)
    imin=fix(0.5*phimin/pi*Np)+1;   % fix rounds towards zero
    imax=fix(0.5*phimax/pi*Np);
    if(phimin < 0)
        imin=1;
    end
    X2=X(imin:imax,:);
    Y2=Y(imin:imax,:);
    Z2=Z(imin:imax,:);
    if(phimin < 0)
        phimin=phimin+2*pi;
        imin=fix(0.5*phimin/pi*Np)+1;
        X=[X(imin:Np-1,:)-2*pi;X2];
        Y=[Y(imin:Np-1,:);Y2];
        Z=[Z(imin:Np-1,:);Z2];
    else
        X=X2;
        Y=Y2;
        Z=Z2;
    end
end

if(Laminar==0 && Machine==1)
    X=(2*pi-X)/pi*180;
end

%------------------------------------------------------------------
% --- Plot --------------------------------------------------------
% Figure
figure;
clf;
axes('FontName',schrift,'FontSize',schrift_size);
if(TypeOfPlot==1)
    pcolor(X,Y,Z);
    %shading flat;
    shading interp;
    if(autoscale==0 || autoscale==2), caxis([min(b) max(b)]); end
else contourf(X,Y,Z,b);
end
set(gca,'layer','top');
set(gca,'linewidth',1.5,'FontName',schrift);
h=findobj('Type','patch'); %nur für contourf
set(h,'EdgeColor','none');
% contourcmap(b,'jet','colorbar','on', ...
%     'FontName',schrift,'FontSize',schrift_size);

% --- additional plots --------------------------------------------
hold on;
% read additional files
for i=plot_manifold:-1:1
    name=char(filenames(i));
    outfile=fopen(name);
    C2=textscan(outfile,'%f%f%*[^\n]','commentStyle','#');
    fclose(outfile);
    phiout=C2{1};
    tout=C2{2};
    if(Laminar==0 && Machine==1), phiout = (2*pi-phiout)/pi*180; end
    plot(phiout,tout,char(styles(i)),'LineWidth',2.5,'markerSize',groesse(i));
end
if(Laminar==1)
    wallfile=fopen('/Users/wingen/c++/d3d/wall.dat');
    C3=textscan(wallfile,'%f%f%f%f','commentStyle','#');
    fclose(wallfile);
    phiout=C3{3};
    tout=C3{4};
    plot(phiout,tout,'k--','LineWidth',2);
end
if(plot_target_limit==1 && Physical ~= 1), plot([min(X(:,1)) max(X(:,1))],[0 0],'k--','LineWidth',1.5); end
if(plot_target_limit==1 && Physical == 1), plot([min(X(:,1)) max(X(:,1))],[-1.228840 -1.228840],'k--','LineWidth',1.5); end
hold off;

% --- Axes limits and labels --------------------------------------
axis tight;
xlim([min(X(:,1)) max(X(:,1))]);
if(Laminar==0 && Target==3 && CameraView==1), ylim([1.371 1.649]);
else ylim([min(Y(1,:)) max(Y(1,:))]);
end
if (Laminar==1)
    xlabel({'R [m]'},'FontName',schrift,'FontSize',schrift_size);
    ylabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
else
    if(Machine==1), xlabel({'$$\phi$$ [deg]'},'Interpreter','latex','FontName',schrift,'FontSize',schrift_size);
    else xlabel({'\fontname{Symbol}j  \fontname{Times}[rad]'},'FontName',schrift,'FontSize',schrift_size);
    end
    if(Physical>=1) 
        if(Target==1) 
            if(Physical==1), ylabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
            else ylabel({'t [cm]'},'FontName',schrift,'FontSize',schrift_size);
            end
        elseif(Target==0), ylabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
        else ylabel({'R [m]'},'FontName',schrift,'FontSize',schrift_size);
        end
    else
        if(Target==1), ylabel({'t_{12}'},'FontName',schrift,'FontSize',schrift_size);
        else ylabel({'t_{23}'},'FontName',schrift,'FontSize',schrift_size);
        end
    end
end

% --- Colormap -----------------------------------------------------
if(autoscale==1), colormap(jet);
else colormap(jet(length(b)-1));
end
if(autoscale==2) 
    MyColorMap=get(gcf,'Colormap');
    if(WhatShallIPlot==2) 
        %MyColorMap(size(MyColorMap,1),:) = 1;
        %MyColorMap(size(MyColorMap,1)-1,:) = 1;
    else MyColorMap(1,:) = 1;
    end
    set(gcf,'Colormap',MyColorMap)
end

% --- set Colorbar, Labels and Ticks ------------------------------
cbar = colorbar('FontName',schrift,'FontSize',schrift_size);
yTickLabelCB = get(cbar, 'YTickLabel');
yTickCB = get(cbar, 'YTick');
%yTickLabelCB = yTickLabelCB(2:end,:);
%yTickCB = yTickCB(2:end);
% cbpos = get(cbar,'Position');

if(WhatShallIPlot==2)
    %text(cbpos(1)+0.5,0.6,'\psi_{Min}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'normalized');
    set(get(cbar,'YLabel'),'String','\fontname{Symbol}y _{\fontname{Times}Min}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90,'VerticalAlignment','bottom')
    %set(cbar,'YTick', [0.901 yTickCB] );
    set(cbar, 'YTickLabel', {yTickLabelCB(1:end-1,:);'SOL'});
elseif(WhatShallIPlot==1)
    %text(cbpos(1)+0.5,0.6,'L_{c} [km]','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'normalized');
    set(get(cbar,'YLabel'),'String','L_{c} [km]','FontName',schrift,'FontSize',schrift_size,'Rotation',-90,'VerticalAlignment','bottom')
    set(cbar,'YTick', [min(b)+0.2*(b(2)-b(1)) yTickCB] );
    set(cbar, 'YTickLabel', {'SOL';yTickLabelCB;});
else
    set(get(cbar,'YLabel'),'String','n_{tor}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90,'VerticalAlignment','bottom')
    %text(cbpos(1)+0.5,0.6,'n_{tor}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'centimeters');
end

%------------------------------------------------------------------
%--- NO export to jpg -----------------------------------
if(printme==0)
% --- Set figure appearance ----------------------------------------
    if(Laminar==0 && Target==1 && Physical~=1), set(gca,'YDir','reverse'); end
    if(Laminar==1), set(gca,'DataAspectRatio',[1 1 1]); end
    if(Laminar==0 && Target==3 && CameraView==1), set(gca,'DataAspectRatio',[42 1 1]); end

%------------------------------------------------------------------
%--- reshape figure and export to jpg -----------------------------------
else
% --- Set figure appearance ----------------------------------------
    if(Laminar==1)
        set(gca, 'Units', 'centimeters');
        set(gca,'DataAspectRatio',[1 1 1]); 
        gcaPBAR = get(gca,'PlotBoxAspectRatio');
        hoehe = 13.5*gcaPBAR(2)/gcaPBAR(1);
        set(gca,'OuterPosition',[0 0 18 hoehe]);

        % set screen size and position of figure window
        set(gcf, 'Units', 'centimeters');
        set(0, 'Units', 'centimeters');
        scrsz = get(0,'ScreenSize');
        set(0, 'Units', 'pixels');
        set(gcf, 'Position', [(scrsz(3)-18)/2 (scrsz(4)-15)/2 18 hoehe]);

        % set size of exported figure
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0.0 0.0 18 hoehe]);

        %cbLabelPosY = 0.5*hoehe;
    else 
        if(Target==1 && Physical~=1), set(gca,'YDir','reverse'); end
        set(gca, 'Units', 'centimeters');
        set(gca, 'OuterPosition', [0 0.5 18 12]);

        % set screen size and position of figure window
        set(gcf, 'Units', 'centimeters');
        set(0, 'Units', 'centimeters');
        scrsz = get(0,'ScreenSize');
        set(0, 'Units', 'pixels');
        set(gcf, 'Position', [(scrsz(3)-18)/2 (scrsz(4)-12)/2 18 12]);

        % set size of exported figure
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0.0 0.0 18 12]);

        if(Machine==1 && VaryPhiRange==0), set(gca, 'XTick', 0:50:350); end
        if(Machine==0 && VaryPhiRange==0)
            set(gca, 'XTick', 0:pi/4:2*pi);
            set(gca, 'XTickLabel', {'0' '' 'pi/2' '' 'pi' '' '3/2 pi' '' '2pi'});
        end

        %cbLabelPosY = 6;
    end
    % Special case
    %     if(Laminar==0 && Target==3 && CameraView==1) 
    %         set(gca,'DataAspectRatio',[42 1 1]); 
    %     end

% --- set Colorbar, Labels and Ticks ------------------------------
%     cbar = colorbar('FontName',schrift,'FontSize',schrift_size);
%     yTickLabelCB = get(cbar, 'YTickLabel');
%     yTickCB = get(cbar, 'YTick');

    set(gca, 'Units', 'centimeters');
    gcapos = get(gca,'Position');
%     if(Laminar==0) cbLabelPosX = gcapos(3) + 3.2;
%     else cbLabelPosX = gcapos(1) + 11.6;
%     end
% 
%     if(WhatShallIPlot==2)
%         text(cbLabelPosX,cbLabelPosY,'\psi_{Min}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'centimeters');
%         set(cbar,'YTick', [0.901 yTickCB 0.999] );
%         set(cbar, 'YTickLabel', {'0.9';yTickLabelCB;'SOL'});
%     elseif(WhatShallIPlot==1)
%         text(cbLabelPosX,cbLabelPosY,'L_{c} [km]','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'centimeters');
%         set(cbar,'YTick', [min(b)+0.6*(b(2)-b(1)) yTickCB 0.395] );
%         set(cbar, 'YTickLabel', {'SOL';yTickLabelCB;'0.4'});
%     else
%         text(cbLabelPosX,cbLabelPosY,'n_{tor}','FontName',schrift,'FontSize',schrift_size,'Rotation',-90, 'Units', 'centimeters');
%     end

% --- Set Colorbar size and Position, reset axes position to match ---
    set(cbar, 'Units', 'centimeters');
    cbpos = get(cbar,'Position');

    if(Laminar==1)
        set(cbar,'location','manual','Position', [11.5+gcapos(1)+correct_pos, cbpos(2), 0.5, cbpos(4)]);
        set(gca, 'Position', [gcapos(1), gcapos(2), gcapos(3), gcapos(4)]);
    else
        set(cbar,'location','manual','Position', [cbpos(1)+0.1, gcapos(2), 0.5, gcapos(4)]);
        if(Machine==1), gcaoff3 = 0.6;
        else gcaoff3 = 0.6;
        end
        set(gca, 'Position', [gcapos(1) gcapos(2) gcapos(3)+gcaoff3 gcapos(4)]);
    end

% --- eport figure to jpg ------------------------------------------
    print('-r300', '-djpeg', filenameout)
    disp(['Figure written to File: ' filenameout '.jpg'])
end

