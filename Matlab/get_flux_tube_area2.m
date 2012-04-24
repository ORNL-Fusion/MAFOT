% clc
% Zuerst d3dfoot.m ausführen!
% Daten aus Workspace werden benutzt
% -- Parameter setzen ----------------------------------------------

% Referenz Punkt
x = 1.545;
y = -1.056;

% Connection length Intervall: [Lmin Lmax]
Lmin = 0.085;
Lmax = 0.125;

% Current density [A/cm^2]
J = 125;

% Zeige Kontroll Plot
show_me = 1;    % 0: Nein   1: Ja

%-- Code -----------------------------------------------------------
Z2 = Z;

%figure(1);

% Alle Flux Tubes
Z2(Z2 > Lmax | Z2 < Lmin) = 0;
[a,b] = find(Z2 > 0);

% Index Position des Referenz-Punktes
[~,xindex]=min(abs(X(:,1)-x));
[value,yindex]=min(abs(Y(1,:)-y));

startindexa = min(find(a-xindex==0));
startindexb = min(find(b-yindex==0));
posa = a(startindexa);
posb = b(startindexb);

% Kontroll Plot vorher
if(show_me == 1)
    figure;
    pcolor(X,Y,Z2)
    shading flat;
    hold on
    %plot(X(posa,posb),Y(posa,posb),'wx','MarkerSize',10,'Linewidth',2);
    hold off
end

% Gehe vom Ref. Punkt zum Rand 
while(Z2(posa+1,posb)>0), posa = posa + 1; end
Z2(posa,posb) = -1; % Markiere Punkt
posa_alt1 = posa;
posb_alt1 = posb;
posa_alt2 = posa;
posb_alt2 = posb;

% Laufe am Rand entlang
Ztest = zeros(5,1);
index = 2;
while(1)
    % Suche Schritt entlang des Randes
    % Suche dabei zuerst in Richtung aus der man kam
    Ztest(mod(1-index+1,4)+1) = Z2(posa+1,posb);
    Ztest(mod(2-index+1,4)+1) = Z2(posa,posb+1);
    Ztest(mod(3-index+1,4)+1) = Z2(posa-1,posb);
    Ztest(mod(4-index+1,4)+1) = Z2(posa,posb-1);
    Ztest(5) = Ztest(1);
    ok = 0;
    for i = 1:4
        if(ok == 0 && Ztest(i) <= 0 && Ztest(i+1) > 0)
            ok = 1;
            % hier Achtung, daß man auch in die richtige Richtung weiter geht: 
            index = mod(i+index+1,4)+1; % nur so klappts
            Ztest(5) = 0;
            tryit = 0;
        end
        if(ok == 1 && sum(Ztest > 0) > 1) % Mehr als ein Schritt noch möglich
            ok = 2;
            posa_alt3 = posa_alt2;
            posb_alt3 = posb_alt2;
            posa_alt2 = posa_alt1;
            posb_alt2 = posb_alt1;
            posa_alt1 = posa;
            posb_alt1 = posb;
            break;
        end
    end
    
    %keinen möglichen Schritt gefunden
    if(ok == 0) 
        if(tryit == 3) 
            break; % Ende
        elseif(tryit == 2) % probiere vor-vor-letzte Position
            posa = posa_alt3;
            posb = posb_alt3;
        elseif(tryit == 1) % probiere vor-letzte Position
            posa = posa_alt2;
            posb = posb_alt2;
        else % probiere letzte Position, an der noch mehrere Schritte möglich waren
            posa = posa_alt1;
            posb = posb_alt1;
        end
        tryit = tryit + 1;
        continue;
    end
    
    % Gehe weiter
    if(index == 1), posb = posb + 1;
    elseif(index == 2), posa = posa - 1;
    elseif(index == 3), posb = posb - 1;
    elseif(index == 4), posa = posa + 1;
    end
    
    % Markiere Punkt
    Z2(posa,posb) = -1; 
    
    % Kontrolliere Verlauf
%     pcolor(X,Y,Z2)
%     shading flat;
%     hold on
%     plot(X(posa,posb),Y(posa,posb),'wx','MarkerSize',10,'Linewidth',2);
%     hold off
%     drawnow
end

% Suche ob noch Punkte im Inneren übersehen wurden
[a,b] = find(Z2 < 0);
for i=3:length(b)-4
    if(b(i) == b(i+1))
        diff = a(i+1) - a(i);
        for j=1:diff-1; Z2(a(i)+j,b(i)) = -1; end
    else continue;
    end
end

% Setze alle unmarkierten auf 0
Z2(Z2 > 0) = 0;
Z2 = -Z2;

% Kontroll Plot nachher
if(show_me == 1)
    figure;
    pcolor(X,Y,Z2)
    shading flat;
    hold on
    %plot(X(posa,posb),Y(posa,posb),'wx','MarkerSize',10,'Linewidth',2);
    hold off
end

% Anzahl der dA's
N = length(find(Z2 > 0));

% Fläche in cm^2 (Nur für Laminar Plot!)
dR = (max(X(:,1)) - min(X(:,1))) / Np;
dZ = (max(Y(1,:)) - min(Y(1,:))) / Nt;
A = N*dR*dZ*100*100
I = J*A

% Get average psi_min inside flux tube
val=C{5};
PSI = reshape(val,Np,Nt);

%PSI = PSI(find(Z2 > 0));
Psi_min_average = sum(PSI(Z2 > 0))/N
