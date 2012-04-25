clear;
old_path = path;
path(old_path,'./functions');
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Code draws particle or penetration depth (and experimental heat flux)
% profile
%
% if(printme == 1) Automatic export to eps
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

FileToOpen = ...%'/home/rack/Masterarbeit/data/jet/079787/foot_out_n1_I36kAt_verylong.dat';
'/.netmount/storage_blast/rack/079459.55918/foot_inn_ef_ft1_9_cam_verylong.dat';
WhatShallIPlot = 2;   % 0: ntor  1: con. length[km]  2: psimin

% experimental data information
FileNameExp = '/home/rack/JET/KL9_data/KL9_T3_79459_55.9214.dat';
x_exp_start = -1.53;
dx_exp = -7/9/100;
N_x_exp = 18;

% calculated data off set
data_offset = -0.1;

% Print option
printme = 1;          % 0: no export to eps file     1: export to eps
praefix = '_nf'; % Comment line if not wanted: Add an additional
                      % string to output filename

% FontName and FontSize
schrift = 'Times';
schrift_size = 20;

% height and width [cm]
fg_height = 12;
fg_width = 16;

% extended options
psi = 271.875;  % position of infra-red camera
deg_avg = 5;    % average over +/- deg_mittel degree
Target = -1;   % 0: Inner  1: Outer  -1: auto mode

% Check for possible file replacement
if(printme==1)
    if(exist('praefix','var')==0), praefix = ''; end
    filenameout = [FileToOpen(1:length(FileToOpen)-4) '_profile' praefix];
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

% automatic mode (read filename and try to identify the data type)
[Laminar, Target] = ident_data_type(FileToOpen, 0, Target);

% read calculated data from file
[X,Y,Z] = read_foot_lam_data(FileToOpen,WhatShallIPlot,0);
X = (2*pi-X)/pi*180;

idx1 = 1;
idx2 = 1;
psi1 = psi-deg_avg;
psi2 = psi+deg_avg;

% find Index
for k=1:length(X(:,1))
    if (abs(X(k,1)-psi1) < abs(X(idx1,1)-psi1))
        idx1 = k;
    end
    if (abs(X(k,1)-psi2) < abs(X(idx2,1)-psi2))
        idx2 = k;
    end
end

% calculate profile
x = Y(X==X(idx2,1));
y = zeros(length(x),1);
for k=idx2:idx1
    y = y + Z(X==X(k,1));
end
y = y / (idx1-idx2+1);

if (Target == 0)
    [dummy x] = t2toroidal(x,0);
else
    [x ~] = t2toroidal(x,1);
end

% read experimental data, adjust data
file = fopen(FileNameExp);
C = textscan(file,'%f%f','commentStyle','#');
fclose(file);
y_exp = C{2};
y_exp = y_exp/1E6;
x_exp = x_exp_start+0.5*dx_exp + dx_exp*(0:N_x_exp-1);
    
% adjust data to the same x-range
if (min(x_exp) < min(x))
    y_exp = y_exp(x_exp >= min(x));
    x_exp = x_exp(x_exp >= min(x));
else
    y = y(x >= min(x_exp));
    x = x(x >= min(x_exp));
end
if (max(x_exp) > max(x))
    y_exp = y_exp(x_exp <= max(x));
    x_exp = x_exp(x_exp <= max(x));
else
    y = y(x <= max(x_exp));
    x = x(x <= max(x_exp));
end

% -------------------------------------------------------------------------
% --- Plot ----------------------------------------------------------------
% Figure
h = figure;
clf;

if (exist('FileNameExp','var'))
    [AX,H1,H2] = plotyy(x,y, x_exp,y_exp);
else
    plot(x,y, 'LineWidth',2)
end

% --- Axis ----------------------------------------------------------------
% advanced options (comment if not wanted)
set(AX,'XGrid','on');

set(AX,'XLim', [min(x),max(x)]);
set(AX(1),'YLim', [min(y)+data_offset,max(y)]);

% add additional Ticks on the left axis because of the data_offset
YTick = get(AX(1),'YTick');
YTick = [YTick(1) + 0.5 * data_offset, YTick];
set(AX(1),'YTick',YTick);

% remove Ticks of the left axis on the right side
set(AX(1),'box','off')
set(AX(2),'XaxisLocation','top','XTickLabel','');

set(H1,'LineWidth',2)
set(H2,'LineWidth',2,'LineStyle','--')

set(AX,'layer','top');
set(AX,'linewidth',1.5,'FontName',schrift,'FontSize',schrift_size);

% Y-axis of simulation
if (WhatShallIPlot==1)
    set(get(AX(1),'Ylabel'),'String','L_{c} [km]', ...
        'FontName',schrift,'FontSize',schrift_size);
else
    set(get(AX(1),'Ylabel'),'String','\psi_{Min}', ...
        'FontName',schrift,'FontSize',schrift_size);
    set(AX(1),'YDir','reverse');
end

% Y-axis of experimental data
set(get(AX(2),'Ylabel'),'String','Heat Flux [MW/m^2]', ...
    'FontName',schrift,'FontSize',schrift_size) 

if (Target==0), 
    xlabel({'Z [m]'},'FontName',schrift, 'FontSize',schrift_size);
else
    xlabel({'R [m]'},'FontName',schrift, 'FontSize',schrift_size);
end

% set screen size and position of figure window
scrsz = get(0,'ScreenSize');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [(scrsz(3)-fg_width)/2 (scrsz(4)-fg_height)/2 ...
    fg_width+2 fg_height+2]);

% set size of exported figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0.0 0.0 fg_width+2 fg_height+2]);


% -------------------------------------------------------------------------
% --- export to eps -------------------------------------------------------
if (printme == 1)
    saveas(h, filenameout, 'epsc')
    disp(['Figure written to File: ' filenameout '.eps'])
end

path(old_path);