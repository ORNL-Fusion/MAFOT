% Zuerst iterfoot.m ausführen!
% Daten aus Workspace werden benutzt
% Nur für ITER footprints
clc
% -- Parameter setzen ----------------------------------------------

% Minimal con. length to consider
cmin = 0.075;

% -- Code ----------------------------------------------------------
% Strikepoint and highest Tip position
if(Target == 1)
    SP = max(Y(find(Z > cmin)))
    Tip = min(Y(find(Z > cmin)))
else
    SP = min(Y(find(Z > cmin)))
    Tip = max(Y(find(Z > cmin)))
end

% Target Plate Coordinates
inner_target = [0.403210000E+01 -0.254550000E+01;
                0.413750000E+01 -0.256500000E+01;
                0.423650000E+01 -0.260590000E+01;
                0.432500000E+01 -0.266650000E+01;
                0.439880000E+01 -0.274420000E+01;
                0.445490000E+01 -0.283560000E+01;
                0.449080000E+01 -0.293660000E+01;
                0.450490000E+01 -0.304290000E+01;
                0.449660000E+01 -0.314970000E+01;
                0.446620000E+01 -0.325260000E+01;
                0.439600000E+01 -0.340160000E+01;
                0.436230000E+01 -0.349580000E+01;   % t<0
                0.417870000E+01 -0.388150000E+01;   % t=0 Point
                0.448570000E+01 -0.390300000E+01];  % t>0
            
outer_target = [0.526550000E+01 -0.425680000E+01;   % t<0
                0.555630000E+01 -0.453180000E+01;   % t=0 Point
                0.555170000E+01 -0.405810000E+01;   % t>0
                0.556090000E+01 -0.397790000E+01;
                0.556210000E+01 -0.381810000E+01;
                0.557870000E+01 -0.371580000E+01;
                0.561110000E+01 -0.361740000E+01;
                0.565840000E+01 -0.352530000E+01;
                0.571960000E+01 -0.344170000E+01;
                0.579310000E+01 -0.336870000E+01;
                0.587720000E+01 -0.330810000E+01;
                0.596960000E+01 -0.326140000E+01;
                0.606830000E+01 -0.322980000E+01;
                0.617060000E+01 -0.321390000E+01;
                0.627420000E+01 -0.321430000E+01];
            
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
if(Target == 1) S = -S; end
S = S*100;  % S in cm

dp = (max(X(:,1)) - min(X(:,1))) / Np /180*pi;  % in rad
dt = (max(Y(1,:)) - min(Y(1,:))) / Nt * 0.01;   % in m

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
    R = tp(index-1,1) + x*(tp(index,1)-tp(index-1,1)); % R(t) in m
    
    % Count Area segments for this t
    N = length(find(Z(:,i) > cmin));
    Area = Area + N*R*dp*dt;    % in m^2
end
Area = Area


