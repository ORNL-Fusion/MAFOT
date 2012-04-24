clear

Target = 1;   % 0: CP  1: Inner + CP  2: Outer  3: Shelf
Physical = 2; % 1: y-Axis in R or Z, depending on target    0: t-Axis   2: t-Axis in cm (inner target only)   

Tmax = 2;       % Temperatur im inneren des Plasmas in [keV]


for k = 1:36
    Ekin = 0.25 + (k-1)*0.5    % Energie in [keV]
    FileToOpen1 = ['foot_in_ip_' num2str(Ekin*100) '.dat'];
    %FileToOpen2 = ['foot_in_ic_' num2str(Ekin*100) '.dat'];

    % --- Read data from input file ---------------------------------------
    file = fopen(FileToOpen1);
    C = textscan(file,'%f%f%f%f%f','commentStyle','#');
    fclose(file);
    L1 = C{4};
    val1 = C{5};
    
%     file = fopen(FileToOpen2);
%     C = textscan(file,'%f%f%f%f%f','commentStyle','#');
%     fclose(file);
%     L2 = C{4};
%     val2 = C{5};
    
    % set X and Y grid only once
    if(k == 1)
        phi = C{1};
        t = C{2};
        sizeofvec = length(t);

        % # grid points of phi
        Np = 1;
        for i = 2:sizeofvec
            if(t(i) == t(i-1)) 
                Np = Np + 1;
            else break;
            end
        end
        Nt = sizeofvec/Np;    % # grid points of t

        X = reshape(phi,Np,Nt);
        Y = reshape(t,Np,Nt);
        
        Z0a = zeros(size(X));
        Z0b = zeros(size(X));
    end

    % Set psi data
    Lc = reshape(L1,Np,Nt);
    Z1 = reshape(val1,Np,Nt);

    %--- Tiefe in psi bestimmen -------------------------------------------
    if(Ekin >= Tmax) 
        psimin = 0.4; 
        psimax = 0.89;
    else
        psimin = newton_T_profil(0.97,1.1*Ekin);     % +-10 Prozent 
        psimax = newton_T_profil(0.98,0.9*Ekin);
    end
    
    %--- Strikepoint position ---------------------------------------------
    SP1 = max(t(L1 >= 0.075));
%    SP2 = max(t(find(L2 >= 0.075)));

    %--- Nur Anteile der entsprechenden Tiefe herausfiltern ---------------
    Z0a(:,:) = 0;    %Z0b(:,:) = 0;
    Z0a(Z1 >= psimin & Z1 <= psimax & Y <= SP1) = Ekin*maxwell(Ekin,Tmax); % anteilig gewichten
%    Z0b(find(Z2 >= psimin & Z2 <= psimax & Y <= SP2)) = Ekin*maxwell(Ekin,Tmax); % anteilig gewichten
%    Z0 = (Z0a + Z0b)/2;
    Z0 = Z0a;
    
    if(k == 1), Zall = Z0; Lcall = Lc;
    else Zall = Zall + Z0; Lcall = Lcall + Lc;
    end
end
Zall = Zall/36;
Lcall = Lcall/36;

% transform t and Y to physical values: length in m
if(Physical >= 1)
    if(Target == 1 && Physical == 1)
        tg = t(t >= 0);
        tl = t(t < 0);
        tg = -1.22884 - tg*0.13756;
        tl = -1.22884 - tl*0.193967;
        t = [tl;tg];
    elseif(Target == 1 && Physical == 2), t = t*19.3967;
    elseif(Target == 0), t = -1.22884 - t*0.193967;
    elseif(Target == 2), t = 1.15285 + t*0.21915;
    elseif(Target == 3), t = 1.372 + t*0.21915;
    end
    Y = reshape(t,Np,Nt);
end

% transform x-Axis to deg and reverse!
X=(2*pi-X)/pi*180;

% Save data to file
save heat_Lc_data X Y Zall Lcall;

% --- Plot --------------------------------------------------------
printme = 0;          % 0: no export to jpg file     1: export to jpg

% (printme==1 only) Comment line if not wanted: Add an additional string to output filename
filenameout = 'foot_in_heat';
%praefix = '_a';

% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Check for possible file replacement
if(printme==1)
    if(exist('praefix','var')==0), praefix = ''; end
    filenameout = [filenameout praefix];
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
xlabel({'$$\phi$$ [deg]'},'Interpreter','latex','FontName',schrift,'FontSize',schrift_size);
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
if(Target==1 && Physical~=1), set(gca,'YDir','reverse'); end

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

    set(gca, 'XTick', 0:50:350);

    %cbLabelPosY = 6;

    set(gca, 'Units', 'centimeters');
    gcapos = get(gca,'Position');

% --- Set Colorbar size and Position, reset axes position to match ---
    set(cbar, 'Units', 'centimeters');
    cbpos = get(cbar,'Position');

    set(cbar,'location','manual','Position', [cbpos(1)+0.1, gcapos(2), 0.5, gcapos(4)]);
    gcaoff3 = 0.6;
    set(gca, 'Position', [gcapos(1) gcapos(2) gcapos(3)+gcaoff3 gcapos(4)]);

% --- eport figure to jpg ------------------------------------------
    print('-r300', '-djpeg', filenameout)
    disp(['Figure written to File: ' filenameout '.jpg'])
end



