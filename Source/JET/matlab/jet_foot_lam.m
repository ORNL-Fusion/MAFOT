clear;
old_path = path;
path(old_path,'./functions');
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Code draws color-contour plot of footprint and laminar data
%
% if(printme == 1) Automatic export to eps
%     
% Warning: in laminar plots check if plot_manifold is either =0
% or only valid data-files are used, otherwise fatal error possible
% (reason unknown!)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --- Set parameter here --------------------------------------------------
% Color scaling in b
WhatShallIPlot = 2;   % 0: ntor  1: con. length[km]  2: psimin
TypeOfPlot = 1;   % 0: contourf  1: pcolor (fast)
autoscale = 2;    % 0: use b    1: scale color automatically 
                  % 2: use b and set SOL value in Colormap to white
b = 0.16:0.005:0.65;   % (con. length only)

% Main Filename
FileToOpen = ...%'/home/rack/Masterarbeit/data/jet/079787/foot_out_n2_I48kAt_verylong.dat';
'/.netmount/storage_blast/rack/079459.55918/foot_out_ef_verylong.dat';

% height and width [cm]
fg_height = 10;
fg_width = 16;

% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Print options
printme = 0;          % 0: no export to eps file     1: export to eps
praefix = '_talk'; % Comment line if not wanted: Add an additional 
                         % string to output filename

% extended options
Laminar = -1;          % 0: Footprint  1: Laminar Plot   -1: auto mode

% (Laminar == 0 only) Specify Target and Axis handling 
Target = -1;    % 0: Inner  1: Outer  -1: auto mode
Physical = 1;   % 0: t-Axis   1: y-Axis in physical length unit, depending on target
Machine = 1;    % 0: x-Axis normal 1: x-Axis reversed and in deg (machine coordinates)
CameraView = 1; % 0: no Modification 1: A rectangle markes the camera observable region

% (Laminar == 0 only) Phi angle area 
VaryPhiRange = 0;  % 0: no Variation 1: Variation 
phimin = (360-65)/180*pi;  %-91
phimax = (360-55)/180*pi; %-43

% Plot additional points 0: none, >=1: Number of Files to plot
plot_manifold = 0;
filenames_mani = {'/.netmount/storage_blast/rack/079459.55918/man_unstr1_ef_0.dat' ...
    '/.netmount/storage_blast/rack/079459.55918/man_stl1_ef_0.dat'};

% Plot additional points from dtstructur files (X,Y,Z,R,phi)
plot_fluxtubemarker = 0;     %0: disabled,   1: enabled
filename_struct = '/.netmount/storage_blast/rack/079459.55918/struct_ef_ft1.dat';
angle = 0;

% Set style for additional data (see Help: LineSpec)
styles_mani = {'-k', '--k', 'mx', 'm.', 'wx', 'w.', 'k+'};
groesse_mani = [10, 10, 10, 15, 10, 15, 6];    % markersizes for additional data 

% (Laminar==1 & printme==1 only) Adjust horizontal colorbar position in [cm], if necessary
correct_pos = 0;    % Default: = 0

% automatic mode (read filename and try to identify the data type)
[Laminar, Target] = ident_data_type(FileToOpen,Laminar, Target);

% Wall parameter
WallFile = '/home/rack/JET/wall.dat';
inn_target_start = 97;
inn_target_end = 80;
out_target_start = 97;
out_target_end = 112;

% positions of breaks for inner target in dimension of R [m], Z [m] and t [cm]
if((Laminar==0))
    wallfile=fopen(WallFile);
    Cwall=textscan(wallfile,'%f%f','commentStyle','#');
    fclose(wallfile);
    Rwall=Cwall{1};
    Zwall=Cwall{2};
    if (inn_target_start > inn_target_end)
        Rdivi = Rwall(inn_target_start:-1:inn_target_end);
        Zdivi = Zwall(inn_target_start:-1:inn_target_end);
    else
        Rdivi = Rwall(inn_target_start:inn_target_end);
        Zdivi = Zwall(inn_target_start:inn_target_end);
    end
    ldivi = sqrt((Rdivi(2:end)-Rdivi(1:end-1)).^2 + ...
        (Zdivi(2:end)-Zdivi(1:end-1)).^2)*100;
    N_tbreaks = length(ldivi)-1;
    tbreaks = zeros(1,N_tbreaks);
    for k=1:N_tbreaks, tbreaks(k) = sum(ldivi(1:k)); end
end

% Parameter for physical Mode
if((Laminar==0) && (Target==0) && (Physical==1))
    tbreaks_start = 8; 
    t2Z_start = tbreaks(tbreaks_start);
    t2Z_end = tbreaks(end);
end

% Specify Data Format
yVariedFirst = 0;     % Default: 0, use 1 only if file format varies y values first

% Specify b for other than con. length
if(WhatShallIPlot == 2), b = 0.95:0.001:1;  % [0.9:0.001:1]
elseif(WhatShallIPlot == 0) 
    b = 4:0.5:40;
    if(autoscale == 2), autoscale = 0; end
end

% plot limit line between CP and inner target: 0=no 1=yes   <-- vielleicht interessant
%if(Laminar==0 && Target==1) plot_target_limit = 1;
%else plot_target_limit = 0;
%end

% Check for possible file replacement
if(printme==1)
    if(exist('praefix','var')==0), praefix = ''; end
    filenameout = [FileToOpen(1:length(FileToOpen)-4) praefix];
    if(exist([filenameout '.eps'],'file')==2)
        disp('File already exists. Enter additional string for Filename.')
        disp('Press just Enter to replace File');
        reply = input('String: ','s');
        if isempty(reply)
            reply = '';
        else filenameout = [filenameout '_' reply];
        end
    end
end

[X,Y,Z] = read_foot_lam_data(FileToOpen,WhatShallIPlot,yVariedFirst);

% transform t to physical values: length in m
sizeofvec = numel(X);
if(Laminar==0 && Physical==1)
    if(Target==1) 
        % only for the first edge!!!
        [rval_used cval_used] = find((Y>=0) & (Y<24.8652));
        if(length(rval_used) ~= sizeofvec) 
            disp('In physical Mode only 0 <= t <= 24.8652 will be plotted!'); 
            rmin = min(rval_used);
            rmax = max(rval_used);
            cmin = min(cval_used);
            cmax = max(cval_used);
            X = X(rmin:rmax,cmin:cmax);
            Y = Y(rmin:rmax,cmin:cmax);
            Z = Z(rmin:rmax,cmin:cmax);
        end
        [Y, ~] = t2toroidal(Y,Target);   %given as R [m]
    else
        [rval_used cval_used] = find((Y>=t2Z_start) & (Y<t2Z_end));
        if(length(rval_used) ~= sizeofvec) 
            disp(['In physical Mode only ' num2str(t2Z_start) ' <= t < ' ...
                num2str(t2Z_end) ' will be plotted!']); 
            rmin = min(rval_used);
            rmax = max(rval_used);
            cmin = min(cval_used);
            cmax = max(cval_used);
            X = X(rmin:rmax,cmin:cmax);
            Y = Y(rmin:rmax,cmin:cmax);
            Z = Z(rmin:rmax,cmin:cmax);
        end
        Y2 = zeros(size(Y));
        for k=tbreaks_start:N_tbreaks-1
            val_used = find((Y>=tbreaks(k)) & (Y<tbreaks(k+1)));
            if (val_used~=0)
                [~, Y2(val_used)] = t2toroidal(Y(val_used),Target);
            end
        end
        Y = reshape(Y2,size(X));     %given as Z [m]
    end
end

% --- Change Z-order at inner target if cameraview is choosen -------------
if ((Laminar==0) && (Physical==1) && (Target==0) && (CameraView==1))
    Y = -Y;
end

% --- change Phi range ----------------------------------------------------
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
    X = (2*pi-X)/pi*180;
end

%--------------------------------------------------------------------------
% --- Plot ----------------------------------------------------------------
% Figure
fg_h = figure;
clf;
if(TypeOfPlot==1)
    pcolor(X,Y,Z);
    %shading flat;
    shading interp;
    if(autoscale==0 || autoscale==2), caxis([min(b) max(b)]); end
else contourf(X,Y,Z,b);
end
set(gca, 'layer','top', 'linewidth',1.5, 'FontName',schrift, ...
    'FontSize',schrift_size);
h = findobj('Type','patch'); %nur für contourf
set(h, 'EdgeColor','none');
% contourcmap(b,'jet','colorbar','on', ...
%     'FontName',schrift,'FontSize',schrift_size);

% --- additional plots ----------------------------------------------------
hold on;
% read additional files
for i=plot_manifold:-1:1
    name = char(filenames_mani(i));
    outfile = fopen(name);
    C2 = textscan(outfile,'%f%f%f%f%f%*[^\n]','commentStyle','#');
    fclose(outfile);
    Rout = C2{4};
    Zout = C2{5};
    plot(Rout,Zout, char(styles_mani(i)), 'LineWidth',2.5, ...
        'markerSize',groesse_mani(i));
end
if (Laminar == 1)
    wallfile = fopen(WallFile);
    C3 = textscan(wallfile,'%f%f','commentStyle','#');
    fclose(wallfile);
    phiout = C3{1};
    tout = C3{2};
    plot(phiout,tout,'k--','LineWidth',2);
else
    if(Target==0)
        if(Physical==0)
            breakmarker = zeros(1,N_tbreaks);
            for k=1:N_tbreaks
                breakmarker(k) = line([X(1) X(end)], ...
                    [tbreaks(k) tbreaks(k)],'LineStyle',':');
            end
            set(breakmarker,'LineWidth',0.5,'Color','k')
        else
            breakmarker = zeros(1,N_tbreaks-tbreaks_start-1);
            for k=tbreaks_start+1:N_tbreaks
                breakmarker(k-tbreaks_start) = line([X(1) X(end)], ...
                    [Zdivi(k) Zdivi(k)],'LineStyle',':');
            end
            set(breakmarker,'LineWidth',0.5,'Color','k')
        end
    end
end
%if(plot_target_limit==1) plot([min(X(:,1)) max(X(:,1))],[0 0],'k--','LineWidth',1.5); end
if (plot_fluxtubemarker == 1)
    file = fopen(filename_struct);
    C_struct = textscan(file,'%f%f%f%f%f','commentStyle','#');
    fclose(file);
    Z_struct = C_struct{3};
    R_struct = C_struct{4};
    phi_struct = C_struct{5};
    plot_struct = find(mod(phi_struct,2*pi) == angle*pi/180);
    plot(R_struct(plot_struct),Z_struct(plot_struct), 'k+', ...
        'LineWidth',2.5, 'markerSize',10);
end
hold off;

% --- Grid, disable if not wanted -----------------------------------------
%grid on;

% --- Axes limits and labels ----------------------------------------------
xlim([min(X(:,1)) max(X(:,1))]);
ylim([min(Y(1,:)) max(Y(1,:))]);
if (Laminar == 1)
    xlabel({'R [m]'}, 'FontName',schrift, 'FontSize',schrift_size);
    ylabel({'Z [m]'}, 'FontName',schrift, 'FontSize',schrift_size);
else
    if(Machine == 1)
        xlabel({'\phi [deg]'}, 'FontName',schrift, ...
            'FontSize',schrift_size);
    else
        xlabel({'\fontname{Symbol}j \fontname{Times}[rad]'}, ...
            'FontName',schrift, 'FontSize',schrift_size);
    end
    if(Physical == 1) 
        if(Target == 0)
            if (CameraView == 1)
                ylabel({'- Z [m]'}, 'FontName',schrift, ...
                    'FontSize',schrift_size);
            else
                ylabel({'Z [m]'}, 'FontName',schrift, ...
                    'FontSize',schrift_size);
            end
        else ylabel({'R [m]'}, 'FontName',schrift, ...
                'FontSize',schrift_size);
        end
    else
        if(Target == 0)
            ylabel({'t_{inn} [cm]'}, 'FontName',schrift, ...
                'FontSize',schrift_size);
        else
            ylabel({'t_{out} [cm]'}, 'FontName',schrift, ...
                'FontSize',schrift_size);
        end
    end
end

% --- Colormap ------------------------------------------------------------
if(autoscale == 1), colormap(jet);
else colormap(jet(length(b)-1));
end
if(autoscale == 2) 
    MyColorMap = get(gcf,'Colormap');
    if(WhatShallIPlot == 2) 
        MyColorMap(size(MyColorMap,1),:) = 1;
        %MyColorMap(size(MyColorMap,1)-1,:) = 1;
        %MyColorMap(size(MyColorMap,1)-2,:) = 1;
        %MyColorMap(size(MyColorMap,1)-3,:) = 1;
        %MyColorMap(size(MyColorMap,1)-4,:) = 1;
        %MyColorMap(size(MyColorMap,1)-5,:) = 1;
        %MyColorMap(size(MyColorMap,1)-6,:) = 1;
    else MyColorMap(1,:) = 1;
    end
    set(gcf,'Colormap',MyColorMap)
end

% --- set Colorbar, Labels and Ticks --------------------------------------
cbar = colorbar('FontName',schrift,'FontSize',schrift_size);
yTickLabelCB = get(cbar, 'YTickLabel');
yTickCB = get(cbar, 'YTick');

if (WhatShallIPlot == 2)
    yTickLabelCB = yTickLabelCB(2:end-1,:);
    yTickCB = yTickCB(2:end-1);
    set(get(cbar,'YLabel'), 'String','\psi_{Min}', 'FontName',schrift, ...
        'FontSize',schrift_size, 'Rotation',-90, ...
        'VerticalAlignment','bottom')
    set(cbar,'YTick', [0.901 yTickCB 0.996] );
    set(cbar, 'YTickLabel', {'0.9';yTickLabelCB;'SOL'});
elseif (WhatShallIPlot == 1)
    %yTickLabelCB = yTickLabelCB(2:end,:);      %remove Tick next to "SOL"
    %yTickCB = yTickCB(2:end);                  %remove Tick next to "SOL"
    set(get(cbar,'YLabel'), 'String','L_{c} [km]', 'FontName',schrift, ...
        'FontSize',schrift_size, 'Rotation',-90, ...
        'VerticalAlignment','bottom')
    set(cbar,'YTick', [min(b)+0.6*(b(2)-b(1)) yTickCB]);
    set(cbar, 'YTickLabel', {'SOL';yTickLabelCB});
else
    set(get(cbar,'YLabel'), 'String','n_{tor}', 'FontName',schrift, ...
        'FontSize',schrift_size, 'Rotation',-90, ...
        'VerticalAlignment','bottom')
end

% --- Draw additional rectangle that denote the cameraview ----------------
if ((Laminar == 0) && (CameraView == 1))
    cv_width = 6;
    cv_pos = 271.875;
    cv_height = abs(Y(1,end)-Y(1,1));
    if ((Physical==1) && (Target==0))
        cv_start = Y(1,end);
    else
        cv_start = Y(1,1);
    end
    line(cv_pos+0.5*cv_width*[-1,1],Y(1,end)*[1,1], 'LineWidth',4, ...
        'Color', 'k')
    rectangle('Position', ...
        [cv_pos-0.5*cv_width, cv_start, cv_width, cv_height], ...
        'LineWidth',2, 'Clipping','off');
    line(cv_pos+0.5*cv_width*[-1,1],Y(1,1)*[1,1], 'LineWidth',4, ...
        'Color', 'k')
end

%--------------------------------------------------------------------------
% --- Set figure appearance -----------------------------------------------
set(gca, 'Units','centimeters');
set(gcf, 'Units','centimeters');
set(0, 'Units','centimeters');
scrsz = get(0, 'ScreenSize');
set(0, 'Units','pixels');
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits','centimeters');

if(Laminar == 1)
    set(gca, 'DataAspectRatio',[1 1 1]); 
    gcaPBAR = get(gca, 'PlotBoxAspectRatio');
    fg_height = fg_width*gcaPBAR(2)/gcaPBAR(1);
    
    %cbLabelPosY = 0.5*hoehe;
else 
    %if(Target==0) set(gca,'YDir','reverse'); end
    
    if(Machine == 1 && VaryPhiRange == 0), set(gca, 'XTick', 0:50:350); end
    if(Machine == 0 && VaryPhiRange == 0)
        set(gca, 'XTick', 0:pi/4:2*pi);
        set(gca, 'XTickLabel', {'0' '' 'pi/2' '' 'pi' '' '3/2 pi' '' ...
            '2pi'});
    end

    %cbLabelPosY = 6;
end

set(gca, 'OuterPosition',[0 0 fg_width fg_height]);

% set screen size and position of figure window
set(gcf, 'Position',[(scrsz(3)-fg_width)/2 (scrsz(4)-fg_height)/2 ...
     fg_width+1 fg_height]);
 
% set size of exported figure
set(gcf, 'PaperPosition', [0.0 0.0 fg_width+1 fg_height]);

% --- Set Colorbar size and Position, reset axes position to match --------
gcapos = get(gca,'Position');
set(cbar, 'Units','centimeters');
cbpos = get(cbar, 'Position');

if (Laminar == 1)
    set(cbar, 'location','manual', 'Position', ...
        [gcapos(3)+0.5+gcapos(1)+correct_pos, cbpos(2), 0.5, cbpos(4)]);
    set(gca, 'Position',[gcapos(1), gcapos(2), gcapos(3), gcapos(4)]);
else
    set(cbar, 'location','manual', ...
        'Position',[cbpos(1)+0.5, gcapos(2), 0.5, gcapos(4)]);
    if(Machine==1), gcaoff3 = 0.6;
    else gcaoff3 = 0.9;
    end
    set(gca, 'Position',[gcapos(1) gcapos(2) gcapos(3)+gcaoff3 gcapos(4)]);
end

% -------------------------------------------------------------------------
% --- export figure to eps ------------------------------------------------
if (printme == 1)
    saveas(fg_h, filenameout, 'epsc')
    disp(['Figure written to File: ' filenameout '.eps'])
end

path(old_path);