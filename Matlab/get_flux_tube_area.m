% clc
% Zuerst d3dfoot.m ausführen!
% Daten aus Workspace werden benutzt
% -- Parameter setzen ----------------------------------------------

% Box innerhalb der gezählt wird
Xmin = 1.512;
Xmax = 1.633;
Ymin = -1.061;
Ymax = -1.046;

% Connection length Intervall: [Lmin Lmax]
Lmin = 0.043;
Lmax = 0.078;

% Current density [A/cm^2]
J = 125;

% Setze Werte ober- bzw. unterhalb einer oder mehrerer Geraden von P1 nach P2 auf Null;
% Länge von use_line bestimmt die Anzahl der Geraden. 
% Punkte im Format P1=[x1,y1,x2,y2,x3,y3,...]
use_line = [0];   % 0: keine Gerade verwenden     1: oberhalb null    -1: unterhalb null
P1 = [1.344 , -1.074 , 1.357 , -1.065, 1.344 , -1.074 , 1.373 , -1.121];
P2 = [1.357 , -1.065 , 1.36 , -1.116 , 1.369 , -1.114 , 1.397 , -1.111];

% Zeige plot nach Anwendung der Geraden (Kontrolle)
show_me = 0;    % 0: Nein   1: Ja

%-- Code -----------------------------------------------------------
Z2 = Z;
for i=1:length(use_line)
    if (use_line(i) == -1)
        a = (P2(2*i) - P1(2*i)) / (P2(2*i-1) - P1(2*i-1));
        b = P1(2*i) - a*P1(2*i-1);
        Z2(find(Y <= a*X + b)) = 0;
    elseif (use_line(i) == 1)
        a = (P2(2*i) - P1(2*i)) / (P2(2*i-1) - P1(2*i-1));
        b = P1(2*i) - a*P1(2*i-1);
        Z2(find(Y >= a*X + b)) = 0;
    end
end

if(show_me == 1)
    figure;
    pcolor(X,Y,Z2)
    shading flat;
end

Z2 = Z2(find(X > Xmin & X < Xmax & Y > Ymin & Y < Ymax));

% Anzahl der dA's
N = length(find(Z2 >= Lmin & Z2 <= Lmax));

% Fläche in cm^2 (Nur für Laminar Plot!)
dR = (max(X(:,1)) - min(X(:,1))) / Np;
dZ = (max(Y(1,:)) - min(Y(1,:))) / Nt;
A = N*dR*dZ*100*100
I = J*A

% Get average psi_min inside flux tube
val=C{5};
PSI=val(1:Np);
for i=2:Nt PSI=[PSI,val((i-1)*Np+1:i*Np)]; end

PSI = PSI(find(X > Xmin & X < Xmax & Y > Ymin & Y < Ymax));
Psi_min_average = sum(PSI(find(Z2 >= Lmin & Z2 <= Lmax)))/N
