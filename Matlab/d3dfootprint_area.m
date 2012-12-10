% Zuerst d3dfoot.m ausführen!
% Daten aus Workspace werden benutzt
% Nur für D3D footprints
clc
% -- Parameter setzen ----------------------------------------------

% Minimal con. length to consider
cmin = 0.075;

% -- Code ----------------------------------------------------------
% Strikepoint and highest Tip position
if(Target == 1)
    SP = max(Y(Z > cmin))   % t in cm; for Z in m swap SP and Tip 
    Tip = min(Y(Z > cmin))  % t in cm; for Z in m swap SP and Tip 
    Extent = abs(SP - Tip)  % in cm
else
    SP = min(Y(Z > cmin))
    Tip = max(Y(Z > cmin))
    Extent = abs(SP - Tip)  % in cm
end

% Target Plate Coordinates in m
inner_target = [1.016 -1.0271;  % t<0
                1.016 -1.223;   % t=0 Point
                1.153 -1.363];  % t>0
            
outer_target = [1.153 -1.363;   % t=0 Point
                1.372 -1.363];  % t>0

            
if(Target == 1) 
    tp = inner_target(end:-1:1,:);
else
    tp = outer_target;
end

% Length along Target Plate
S = zeros(length(tp),1);
S(1) = -norm(tp(2,:)-tp(1,:));
S(2) = 0;
for i=3:length(tp)
    S(i) = S(i-1) + norm(tp(i,:)-tp(i-1,:));
end
if(Target == 1), S = -S; end
S = S*100;  % S in cm

dp = (max(X(:,1)) - min(X(:,1))) / Np /180*pi;  % in rad
dt = (max(Y(1,:)) - min(Y(1,:))) / Nt;          % in cm

% Calculate Area of Footprint
Area = 0;
for i = 1:Nt
    t = Y(1,i); % in cm

    % Calculate R(t)
    if(Target == 1)
        index = min(find(S < t));   % t between S(index-1) and S(index)
    else
        index = min(find(S > t));
    end
    x = (S(index-1) - t)/(S(index-1) - S(index));
    R = (tp(index-1,1) + x*(tp(index,1)-tp(index-1,1)))*100; % R(t) in cm
    
    % Count Area segments for this t
    N = length(find(Z(:,i) > cmin));
    Area = Area + N*R*dp*dt;    % in cm^2
end
Area = Area


