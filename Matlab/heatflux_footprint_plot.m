% Run heatflux_footprint_generate.m first!
%--------------------------------------------------------------------------
% or load data from file
load heat_Lc_data
Target = 1;   % 0: CP  1: Inner + CP  2: Outer  3: Shelf
Physical = 2; % 1: y-Axis in R or Z, depending on target    0: t-Axis   2: t-Axis in cm (inner target only)   
%--------------------------------------------------------------------------

printme = 0;          % 0: no export to jpg file     1: export to jpg

% (printme==1 only) Comment line if not wanted: Add an additional string to output filename
filenameout = 'foot_in_heat';
%praefix = '_a';

% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Check for possible file replacement
if(printme==1)
    if(exist('praefix')==0) praefix = ''; end
    filenameout = [filenameout praefix];
    if(exist([filenameout '.jpg'])==2)
        disp('File already exists. Enter additional string for Filename.')
        disp('Press just Enter to replace File');
        reply = input('String: ','s');
        if isempty(reply)
            reply = '';
        else filenameout = [filenameout '_' reply];
        end
    end
end

% Figure
figure;
clf;
axes('FontName',schrift,'FontSize',schrift_size);
pcolor(X,Y,Zall);
shading interp;
set(gca,'layer','top');
set(gca,'linewidth',1.5,'FontName',schrift);

% --- Axes limits and labels --------------------------------------
axis tight;
xlim([min(X(:,1)) max(X(:,1))]);
ylim([min(Y(1,:)) max(Y(1,:))]);
xlabel({'\phi [deg]'},'FontName',schrift,'FontSize',schrift_size);
if(Physical>=1)
    if(Target==1)
        if(Physical==1) ylabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
        else ylabel({'t [cm]'},'FontName',schrift,'FontSize',schrift_size);
        end
    elseif(Target==0) ylabel({'Z [m]'},'FontName',schrift,'FontSize',schrift_size);
    else ylabel({'R [m]'},'FontName',schrift,'FontSize',schrift_size);
    end
else
    if(Target==1) ylabel({'t_{12}'},'FontName',schrift,'FontSize',schrift_size);
    else ylabel({'t_{23}'},'FontName',schrift,'FontSize',schrift_size);
    end
end
if(Target==1 && Physical~=1) set(gca,'YDir','reverse'); end

% --- Colormap -----------------------------------------------------
colormap(jet);
MyColorMap=get(gcf,'Colormap');
MyColorMap(1,:) = 1;
set(gcf,'Colormap',MyColorMap)

% --- set Colorbar, Labels and Ticks ------------------------------
cbar = colorbar('FontName',schrift,'FontSize',schrift_size);
yTickLabelCB = get(cbar, 'YTickLabel');
yTickCB = get(cbar, 'YTick');
yTickLabelCB = yTickLabelCB(2:end,:);
yTickCB = yTickCB(2:end);

set(get(cbar,'YLabel'),'String','T [keV]','FontName',schrift,'FontSize',schrift_size,'Rotation',-90,'VerticalAlignment','bottom')
% set(cbar,'YTick', [0.901 yTickCB 0.999] );
% set(cbar, 'YTickLabel', {'0.9';yTickLabelCB;'SOL'});

%------------------------------------------------------------------
%--- NO export to jpg -----------------------------------
if(printme==0)
% --- Set figure appearance ----------------------------------------

%------------------------------------------------------------------
%--- reshape figure and export to jpg -----------------------------------
else
% --- Set figure appearance ----------------------------------------
    set(gca, 'Units', 'centimeters');
    set(gca, 'OuterPosition', [0 0 18.5 12.5]);

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

    set(gca, 'XTick', [0:50:350]);

    %cbLabelPosY = 6;

    set(gca, 'Units', 'centimeters');
    gcapos = get(gca,'Position');

% --- Set Colorbar size and Position, reset axes position to match ---
    set(cbar, 'Units', 'centimeters');
    cbpos = get(cbar,'Position');

    set(cbar,'location','manual','Position', [cbpos(1)+0.5, gcapos(2), 0.5, gcapos(4)]);
    gcaoff3 = 0.6;
    set(gca, 'Position', [gcapos(1) gcapos(2) gcapos(3)+gcaoff3 gcapos(4)]);

% --- eport figure to jpg ------------------------------------------
    print('-r300', '-djpeg', filenameout)
    disp(['Figure written to File: ' filenameout '.jpg'])
end


