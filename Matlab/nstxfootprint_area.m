% Zuerst nstxfoot.m ausführen!
% Daten aus Workspace werden benutzt
% Nur für D3D footprints
clc
% -- Parameter setzen ----------------------------------------------

% Minimal con. length to consider
cmin = 0.046;

% -- Code ----------------------------------------------------------
% Strikepoint and highest Tip position
if(Target == 1)
    SP = min(Y(Z > cmin))   % Z in m
    Tip = max(Y(Z > cmin))  % Z in m
    Extent = abs(SP - Tip)*100  % in cm
else
    SP = min(Y(Z > cmin))   % R in m
    Tip = max(Y(Z > cmin))  % R in m
    Extent = abs(SP - Tip)*100  % in cm
end

% Target Plate Coordinates in m; due to symmetry only one set required
% upper inner target plate, length 40.66 cm
%inner_target_up = [0.2794 1.1714;
%                   0.2794 1.578];
% upper outer target plate, length 27.33 cm
%outer_target_up = [0.2979 1.6034;
%                   0.5712 1.6034];
% lower inner target plate, length 40.66 cm
inner_target = [0.2794 -1.1714;
                0.2794 -1.578];
% lower outer target plate, length 27.33 cm
outer_target = [0.2979 -1.6034;
                0.5712 -1.6034];
            
if(Target == 1) 
    tp = inner_target;
else
    tp = outer_target;
end

dp = (max(X(:,1)) - min(X(:,1))) / Np /180*pi;  % in rad
dt = (max(Y(1,:)) - min(Y(1,:))) / Nt * 100;    % in cm

% Calculate Area of Footprint
Area = 0;
for i = 1:Nt
    % Calculate R(i) in cm
    if(Target == 1)
        R = tp(1,1)*100;
    else
        R = Y(1,i)*100;
    end
    
    % Count Area segments for this i
    N = length(find(Z(:,i) > cmin));
    Area = Area + N*R*dp*dt;    % in cm^2
end
Area = Area


