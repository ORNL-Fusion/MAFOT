% Parameter
%--------------
useFiles = 1;           % 1: read from files    0: use workspace
FileToOpen1 = 'foot_inup.dat';
%FileToOpen2 = 'foot_outdwn.dat';   % Comment, if not wanted

useLogScale = 1;        %0: linear y-Axis   1: log(y+1) y-Axis

M = 200;                % Number of bars
SOL_Range_Lc = 0.046;   % Minimum connection length; D3D = 0.075; ITER = 0.22; NSTX = 0.046
SOL_Range_psi = 0.999;  % Maximum penetration depth; D3D = 0.997; ITER = 0.997; NSTX = 0.999
Lc_cutoff = 2.6;        % Cutoff Lc x-Axis at value [km], if value > max(Lx): no cutoff

Lc_ft = 0.085;           % max Flux Tube con. Length = 2pol. Turns; D3D = 0.13; ITER = 0.44; NSTX = 0.085
Lc_long = 0.2;          % another con. Length; D3D = 0.4; ITER = 1; NSTX = 0.2
correct_psi = 1;
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

    if(exist('FileToOpen2','var')==1)
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
    WhatToPlot = 0; % 0: Plot both, 1: Lc, 2: psimin

    % Correct psi: use Lc to set psi = 0 in SOL
    if(correct_psi == 1), Psi(Lc <= SOL_Range_Lc) = 1; end
    
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
        if(k(i) > M), continue; end  % Values larger than zmax are ignored
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
        if(k(i) <= 0), continue; end  % All values smaller than zmin are ignored
        if(k(i) > M), k(i) = M; end   % All values greater than zmax are set to zmax
        L(k(i)) = L(k(i)) + 1;
        count = count + 1;
    end

    % Normalize
    L = L/count*100;

    % Lc Axis
    Lx = zmin + (1:M)'*dz - dz/2;
    
    % Statistical data
    Lc_average = round(sum(Lc(Lc>zmin))/length(Lc(Lc>zmin))*100)/100;
    flf1b = round(10*sum(L(1:max(find(Lx<=Lc_ft)))))/10; 
    flf2b = round(10*sum(L(1:max(find(Lx<=Lc_long)))))/10;
    
    % Cutoff
    L = L(Lx<=Lc_cutoff);
    Lx = Lx(Lx<=Lc_cutoff);
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
            'YTick',[1,1.2,1.4,1.6,1.8,2,3,4,5,6,8,11,16,21,31,41,51,71],...
            'YTickLabel',{'0';'0.2';'0.4';'0.6';'0.8';'1';'2';'3';'4';'5';'7';'10';'15';'20';'30';'40';'50';'70'},...
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

    text(min(Px)+4*dz,0.92*max(P),str,...
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
    str = {['<L_{c}> [km] = ' num2str(Lc_average) ' km']; ...
        ['flf_{' num2str(Lc_ft*1000) 'm} = ' num2str(flf1b) '%']; ...
        ['flf_{' num2str(Lc_long*1000) 'm} = ' num2str(flf2b) '%']};

    text(max(Lx)-4*dz,0.92*max(L),str,...
        'BackgroundColor','w','FontName',schrift,...
        'FontSize',schrift_size,'VerticalAlignment','top',...
        'HorizontalAlignment','right');
end

figure(gcf) % Bring current figure window to the front

