% Zuerst d3dfoot2.m ausführen!
% Daten aus Workspace werden benutzt
% -- Parameter setzen ----------------------------------------------

clc

% Minimal con. length to consider
cmin = 0.22;

% Sample
channels = 100;
laenge = 10;
dS = -laenge / channels;

% Strikepoint position
SP = max(Y(find(Z > cmin)));

area_ratio = zeros(channels,1);

for i = 1:channels
    area_ratio(i) = length(find(Z(Y <= SP+(i-1)*dS & Y > SP+i*dS) > cmin)) / length(find(Y <= SP+(i-1)*dS & Y > SP+i*dS));    
end
area_ratio