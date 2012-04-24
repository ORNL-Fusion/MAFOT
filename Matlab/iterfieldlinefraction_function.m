function [psi_min,flf1a,flf2a,Lc_average,flf1b,flf2b] = iterfieldlinefraction_function(FileToOpen1);
% Parameter
%--------------
useFiles = 1;           % 1: read from files    0: use workspace
%FileToOpen1 = 'foot_out_n3_45kA_86-0-34.dat';
%FileToOpen2 = 'foot_in_n3_45kA_86-0-34.dat';   % Comment, if not wanted

useLogScale = 1;        %0: linear y-Axis   1: log(y+1) y-Axis

M = 200;                % Number of bars
SOL_Range_Lc = 0.22;    % Minimum connection length
SOL_Range_psi = 0.997;  % Maximum penetration depth
Lc_cutoff = 16.5;       % Cutoff Lc x-Axis at value [km], if value > max(Lx): no cutoff

%------------------------------------------------------------------
%------------------------------------------------------------------
% Code draws histogramms for Connection length and Penetration depth
% specifically designed for footprint data
%
% You can eiter read in one or two footprint files
% then you get both histogramms
%
% or use workspace data from d3dfoot.m (ITER or NSTX as well)
% and get the one corresponding to the data in 'Z'
% 'Z' and 'WhatShallIPlot' are used from workspace
%
% useLogScale aktivates a logarithmic scale for better visibility
% to avoid negative numbers, log(data + 1) is ploted
% the axis tick labes are manipulated so that the values are correct
%------------------------------------------------------------------
%------------------------------------------------------------------

% Read Files
%--------------
if(useFiles==1)
    file=fopen(FileToOpen1);
    C=textscan(file,'%f%f%f%f%f','commentStyle','#');
    fclose(file);
    Lc1=C{4};
    Psi1=C{5};

    if(exist('FileToOpen2')==1)
        file=fopen(FileToOpen2);
        C=textscan(file,'%f%f%f%f%f','commentStyle','#');
        fclose(file);
        Lc2=C{4};
        Psi2=C{5};

        Psi = [Psi1;Psi2];
        Lc = [Lc1;Lc2];
    else
        Psi = Psi1;
        Lc = Lc1;
    end
    WhatToPlot = 0; % Plot both
    
% Take Data from Workspace
%--------------
else
    N = size(Z,1)*size(Z,2);
    if(WhatShallIPlot == 2)
        Psi = reshape(Z,N,1);
    else
        Lc = reshape(Z,N,1);
    end
    WhatToPlot = WhatShallIPlot;
end

% Psi Histogramm
%--------------
if(WhatToPlot ~= 1)
    zmin = min(Psi)-0.01;
    zmax = SOL_Range_psi;
    dz = (zmax-zmin)/M;
    N = length(Psi);

    P = zeros(M,1);
    k = ceil((Psi-zmin)/dz);
    count = 0;
    for i = 1:N
        if(k(i) > M) continue; end  % Values larger than zmax are ignored
        if(k(i) <= 0) continue; end  % All values smaller than zmin are ignored (does not occur usually)
        P(k(i)) = P(k(i)) + 1;
        count = count + 1;
    end

    % Normalize
    P = P/count*100;

    % psi Axis
    Px = zmin + (1:M)'*dz - dz/2;
    
    % Statistical data
    psi_min = round(min(Psi)*1000)/1000;
    flf1a = 100-round(10*sum(P(1:max(find(Px<=0.98)))))/10;
    flf2a = 100-round(10*sum(P(1:max(find(Px<=0.95)))))/10;
end

% Lc Histogramm
%--------------
if(WhatToPlot ~= 2)
    zmin = SOL_Range_Lc;
    zmax = max(Lc);
    dz = (zmax-zmin)/M;
    N = length(Lc);

    L = zeros(M,1);
    k = ceil((Lc-zmin)/dz);
    count = 0;
    for i = 1:N
        if(k(i) <= 0) continue; end  % All values smaller than zmin are ignored
        if(k(i) > M) continue; end  % All values larger than zmax are ignored (does not occur usually)
        L(k(i)) = L(k(i)) + 1;
        count = count + 1;
    end

    % Normalize
    L = L/count*100;

    % Lc Axis
    Lx = zmin + (1:M)'*dz - dz/2;
    
    % Statistical data
    Lc_average = round(sum(Lc(find(Lc>zmin)))/length(Lc(find(Lc>zmin)))*100)/100;
    flf1b = round(10*sum(L(1:max(find(Lx<=0.44)))))/10;
    flf2b = round(10*sum(L(1:max(find(Lx<=1)))))/10;
    
    % Cutoff
    L = L(find(Lx<=Lc_cutoff));
    Lx = Lx(find(Lx<=Lc_cutoff));
end

% Output
%--------------
% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Figure Psi
%--------------
if(WhatToPlot ~= 1)
    figure;
    clf
    axes('FontName',schrift,'FontSize',schrift_size);
    
    % Plot
    if(useLogScale == 1)
        P = P+1;
        bar(Px,P,'BaseValue',1);
    else 
        bar(Px,P);
    end
    set(gca,'linewidth',1.5,'FontName',schrift);

    % Axis Limits and Labels
    dz = Px(2)-Px(1);
    xlim([min(Px)-dz 1]);
    if(useLogScale == 1)
        ylim([0 1.1*max(P)]);
        set(gca,'YScale','log',...
            'YTick',[1,2,3,4,5,6,8,11,16,21,31,41,51,71],...
            'YTickLabel',{'0';'1';'2';'3';'4';'5';'7';'10';'15';'20';'30';'40';'50';'70'},...
            'YMinorTick','off');
        set(gca,'XGrid','on','YGrid','on','YMinorGrid','off');
    else 
        ylim([-0.05*max(P) 1.1*max(P)]);
        grid on;
    end
    xlabel({'\fontname{Symbol}y _{\fontname{Times}Min}'},'FontName',schrift,'FontSize',schrift_size);
    ylabel({'field line fraction [%]'},'FontName',schrift,'FontSize',schrift_size);

    % Inlet
    str = {['\fontname{Symbol}y _{\fontname{Times}Min} \fontname{Times} = ' num2str(psi_min)]; ...
        ['flf_{98} = ' num2str(flf1a) '%']; ...
        ['flf_{95} = ' num2str(flf2a) '%']};

    text(min(Px)+2*dz,1.05*max(P),str,...
        'BackgroundColor','w','FontName',schrift,...
        'FontSize',schrift_size,'VerticalAlignment','top');
end

% Figure Lc
%--------------
if(WhatToPlot ~= 2)
    figure;
    clf
    axes('FontName',schrift,'FontSize',schrift_size);

    % Plot
    if(useLogScale == 1)
        L = L+1;
        bar(Lx,L,'BaseValue',1);
    else 
        bar(Lx,L);
    end
    set(gca,'linewidth',1.5,'FontName',schrift);

    % Axis Limits and Labels
    dz = Lx(2)-Lx(1);
    xlim([min(Lx)-dz max(Lx)+dz]);
    if(useLogScale == 1)
        ylim([0 1.1*max(L)]);
        set(gca,'YScale','log',...
            'YTick',[1,2,3,4,5,6,8,11,16,21,31,41,51,71],...
            'YTickLabel',{'0';'1';'2';'3';'4';'5';'7';'10';'15';'20';'30';'40';'50';'70'},...
            'YMinorTick','off');
        set(gca,'XGrid','on','YGrid','on','YMinorGrid','off');
    else 
        ylim([-0.05*max(L) 1.1*max(L)]);
        grid on;
    end
    xlabel({'L_c [km]'},'FontName',schrift,'FontSize',schrift_size);
    ylabel({'field line fraction [%]'},'FontName',schrift,'FontSize',schrift_size);

    % Inlet
    str = {['\langle L_{c}\rangle [km] = ' num2str(Lc_average) ' km']; ...
        ['flf_{440m} = ' num2str(flf1b) '%']; ...
        ['flf_{1km} = ' num2str(flf2b) '%']};

    text(max(Lx)-2*dz,1.05*max(L),str,...
        'BackgroundColor','w','FontName',schrift,...
        'FontSize',schrift_size,'VerticalAlignment','top',...
        'HorizontalAlignment','right');
end



