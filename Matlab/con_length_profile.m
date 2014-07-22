
% Execute d3dfoot.m first!
% Data from Workspace is used
%
% Con. length is averaged over angle width
%
%---------------- Set parameter here -------------------------------

% Angle and Average range of profile in deg
phi = 0;
dphi = 5;      % average over [phi-dphi,phi+dphi]

% Type of Profile
profile = 0;    % 0 = Lc,   1 = HeatFlux

% Strike-Point Position
SP = 0;
shift = 0;      % additional shift of Profile

% Save Data to file
saveme = 1;
shot = '132741';
time = '3950';
if(profile == 1) savename = ['./heatflux_' num2str(phi) 'deg_' shot '_' time '_inner_ip.dat'];
else savename = ['./con_length_' num2str(phi) 'deg_nopr.dat'];
end

%-------------------------------------------------------------------
%-------------------- Code -----------------------------------------
% get Profile
%x = (1.0161+Y(1,:)*(1.15285-1.0161)-0.03)'; %SP - Y(1,:);
if(profile == 1) Z = Zall; end

index = find(X(:,1)>phi-dphi & X(:,1)<phi+dphi);
Zp = Z(min(index):max(index),:);
x = Y(1,:)';
y = (sum(Zp,1)./size(Zp,1))';
%x = SP - x + shift;

% save Variables to disk
if(saveme == 1) 
    A = [x y]; 
    save(savename,'A','-ascii'); 
end

figure
plot(x,y,'k-', 'linewidth',1.5);
xlim([min(x) max(x)]);
xlim([1.4 1.6]);
ylim([0 0.2]);
xlabel({'R [m]'});
ylabel({'L_c [km]'});
grid on
%close all;

figure(gcf) % Bring current figure window to the front
