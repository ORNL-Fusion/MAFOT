% Parameter
%--------------
FileToOpen = 'lam_psi_nopr_4140.dat';

useLogScale = 0;        %0: linear y-Axis   1: log(y+1) y-Axis

M = 100;                % Number of bars
SOL_Range_Lc = 0.075;   % Minimum connection length
SOL_Range_psi = 1.02;   % Maximum normalized flux

%------------------------------------------------------------------
%------------------------------------------------------------------
% Code draws histogramms for Connection length and Penetration depth
% specifically designed for laminar_psi data
%
% You can eiter read in a file an get both histograms
% or use workspace data from d3dfoot.m (ITER or NSTX as well)
% and get the one corresponding to the data in 'Z'
% 'Z' and 'WhatShallIPlot' are used from workspace
%
% useLogScale aktivates a logarithmic scale for better visibility
% to avoid negative numbers, log(data + 1) is ploted
% the axis tick labes are manipulated so that the values are correct
%------------------------------------------------------------------
%------------------------------------------------------------------

% Read File
%--------------
file=fopen(FileToOpen);
C=textscan(file,'%f%f%f%f%f%f%*[^\n]','commentStyle','#');
fclose(file);
Psi = C{2};
N_tor = C{3};
Lc = C{4};
Psimin = C{5};
Psimax = C{6};

Nmax = max(N_tor);
    
% Psi Histogramm
%--------------
zmin = min(Psi) - 0.01;
zmax = SOL_Range_psi;
dz = (zmax-zmin)/M;
N = length(Psi);
%dp = (max(Psi) - min(Psi))/(Np);

P = zeros(M,1);
k = ceil((Psi-zmin)/dz);
for i = 1:N
    if(k(i) > M), continue; end  % Values larger than zmax are ignored
    if(N_tor(i) == Nmax), continue; end % use only lost field lines
    P(k(i)) = P(k(i)) + 1;
end

% Normalize
for m = 1:M
    if(~isempty(k(k == m))), P(m) = P(m)/length(k(k == m))*100; end
end

% psi Axis
Px = zmin + (1:M)'*dz - dz/2;

% Statistical data
psi_min = round(min(Psimin((N_tor < Nmax) & (Psi < 1)))*1000)/1000;

% Output
%--------------
% FontName and FontSize
schrift = 'Times';
schrift_size = 18;

% Figure Psi
%--------------
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
        'YTick',[1,2,3,4,5,6,8,11,16,21,31,41,51,71,101],...
        'YTickLabel',{'0';'1';'2';'3';'4';'5';'7';'10';'15';'20';'30';'40';'50';'70';'100'},...
        'YMinorTick','off');
    set(gca,'XGrid','on','YGrid','on','YMinorGrid','off');
else 
    ylim([-0.05*max(P) 1.1*max(P)]);
    grid on;
end
xlabel({'\fontname{Symbol}y'},'FontName',schrift,'FontSize',schrift_size);
ylabel({'field line loss fraction [%]'},'FontName',schrift,'FontSize',schrift_size);

% Inlet
str = {['\fontname{Symbol}y _{\fontname{Times}Min} \fontname{Times} = ' num2str(psi_min)];};

text(min(Px)+4*dz,0.92*max(P),str,...
    'BackgroundColor','w','FontName',schrift,...
    'FontSize',schrift_size,'VerticalAlignment','top');

figure(gcf) % Bring current figure window to the front

