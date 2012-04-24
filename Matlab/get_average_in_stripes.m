% Zuerst d3dfoot.m ausführen!
% Daten aus Workspace werden benutzt
% -- Parameter setzen ----------------------------------------------

clc
show_figures = 1; % 0: alle Figuren werden am Ende geschlossen   1: Nichts passiert

% Box innerhalb der gezählt wird
xmin = 209;
xmax = 280;

% Connection length bzw. psi Intervall
cmin = 0.0;               % 0.0;
cmax = 0.997;                 % 0.998;    % value <= max

% Anzahl der Stripes die untersucht werden
stripes = 2;

%----------------------------------------------------------
%----------------------------------------------------------
% Code berechnet Mittelwert innerhalb der einzelnen Streifen
% psi_min oder Con. length möglich
%
% Code berechnet außerdem die t-Werte für den oberen und unteren Rand
% jedes Streifens
%
% Benötigt Workspace Daten von d3dfoot.m
%
%----------------------------------------------------------
%----------------------------------------------------------

% Minimiere X,Y,Z auf relevante Box
X2 = X(:,1);
index = find(X2 >= xmin & X2 <= xmax);

X2 = X(index,:);
Y2 = Y(index,:);
Z2 = Z(index,:);
[n,m] = size(Z2);

% Setze Werte außerhalb des Intervalls auf null
Z2(find(Z2 > cmax | Z2 < cmin)) = 0;
figure;
clf;
pcolor(X2,Y2,Z2);
shading flat;

% Finde stripes
ivec1 = find(Z2(1,:) > 0);
ivecn = find(Z2(n,:) > 0);

% Finde y-Coordinate für Beginn und Ende der einzelnen Streifen in der
% Mitte der Box
ivec = find(Z2(round(n/2),:) > 0);
iveca = [ivec(1)+1,ivec];   % erster und letzter wert in ivec werden nicht gefunden
ivecb = [ivec,1];

t_start = [Y(1,ivec(1)),Y(1,ivec(find(ivecb-iveca > 1)))];
t_end = [Y(1,ivec(find(ivecb-iveca > 1)-1)), Y(1,ivec(length(ivec)))];

t_start = t_start(1:stripes)';
t_end = t_end(1:stripes)';
disp('      t_start                 t_end')
disp([t_start    t_end])
disp('Strikeline')
disp(Y(1,ivec(length(ivec))))

% Prüfe ob Stripes über volle Breite, evtl. kompensiere
X_index = n-1;
while(abs(length(ivec1)-length(ivecn)) > 4 & X_index > 2)
    X_index = round(0.5*X_index);
    ivec1 = find(Z2(n-X_index,:) > 0);
end

% Berechne Start und Endpunkte für Geraden
P1 = zeros(2,stripes+1);
j = 1;
for i = 1:length(ivec1)-1
    if(abs(ivec1(i+1)-ivec1(i)) > 3)
        index = round(0.5*(ivec1(i+1)+ivec1(i)));
        P1(1,j) = X2(n-X_index,index);
        P1(2,j) = Y2(n-X_index,index);
        j = j + 1;
        if(j > stripes) break; end;
    end
end

Pn = zeros(2,stripes+1);
j = 1;
for i = 1:length(ivecn)-1
    if(abs(ivecn(i+1)-ivecn(i)) > 3)
        index = round(0.5*(ivecn(i+1)+ivecn(i)));
        Pn(1,j) = X2(n,index);
        Pn(2,j) = Y2(n,index);
        j = j + 1;
        if(j > stripes) break; end;
    end
end

% Letzte Gerade (suche von hinten neu starten)
j = stripes + 1;
for i = length(ivec1):-1:2
    if(abs(ivec1(i-1)-ivec1(i)) > 3)
        index = round(0.5*(ivec1(i-1)+ivec1(i)));
        P1(1,j) = X2(n-X_index,index);
        P1(2,j) = Y2(n-X_index,index);
        break;
    end
end

for i = length(ivecn):-1:2
    if(abs(ivecn(i-1)-ivecn(i)) > 3)
        index = round(0.5*(ivecn(i-1)+ivecn(i)));
        Pn(1,j) = X2(n,index);
        Pn(2,j) = Y2(n,index);
        break;
    end
end
        
% Mittele über Streifen
values = zeros(stripes,1);
Z4 = Z2;
for i = 1:stripes
    Z3 = Z4;
    a = (Pn(2,i) - P1(2,i)) / (Pn(1,i) - P1(1,i));
    b = P1(2,i) - a*P1(1,i);
    Z3(find(Y2 > a*X2 + b)) = 0;    % Setze Werte oberhalb der Geraden auf null

    figure;
    clf;
    pcolor(X2,Y2,Z3);
    shading flat;
    hold on;
    plot(P1(1,i),P1(2,i),'wx','markerSize',8);
    plot(Pn(1,i),Pn(2,i),'wx','markerSize',8);
    hold off;
    
    N = length(find(Z3 > 0))
    Average = sum(sum(Z3))/N
    values(i) = Average;

    Z4(find(Z3 > 0)) = 0;
end

% Bestimme mittlere Eindringtiefe nahe am Strike point
Z3 = Z2;
i = stripes;  % vorletzte Gerade
a = (Pn(2,i) - P1(2,i)) / (Pn(1,i) - P1(1,i));
b = P1(2,i) - a*P1(1,i);
Z3(find(Y2 < a*X2 + b)) = 0;    % Setze Werte unterhalb der Geraden auf null

i = stripes + 1;  % letzte Gerade
a = (Pn(2,i) - P1(2,i)) / (Pn(1,i) - P1(1,i));
b = P1(2,i) - a*P1(1,i);
Z3(find(Y2 > a*X2 + b)) = 0;    % Setze Werte oberhalb der Geraden auf null

figure;
clf;
pcolor(X2,Y2,Z3);
shading flat;
hold on;
plot(P1(1,i),P1(2,i),'wx','markerSize',8);
plot(Pn(1,i),Pn(2,i),'wx','markerSize',8);
hold off;

N_rest = length(find(Z3 > 0))
Average_rest = sum(sum(Z3))/N_rest

if(show_figures == 0) close all; end

