clear;

% Z Skalierung in b
b=[0:1:20];

% Phi Winkel Bereich
phimin=0;
phimax=2*pi;
namen=['foot_2250foot_low.dat';
    'foot_2350foot_low.dat';
    'foot_2450foot_low.dat';
    'foot_2550foot_low.dat';
    'foot_2690foot_low.dat';
    'foot_2790foot_low.dat';
    'foot_2890foot_low.dat';
    'foot_2965foot_low.dat';
    'foot_3050foot_low.dat';
    'foot_3150foot_low.dat';
    'foot_3250foot_low.dat';
    'foot_3450foot_low.dat'];

%mov = avifile('test.avi')
scrsz = get(0,'ScreenSize');
figure('Position',[0.05*scrsz(3),0.1*scrsz(4),0.9*scrsz(3),0.8*scrsz(4)])
for j=1:length(namen(:,1))
% Filename spezifizieren
file=fopen(namen(j,:));

% Einlesen der Daten aus Main File
% die Matrizen hier sind transponiert zu den unteren, ist aber egal
C=textscan(file,'%f%f%f%f','commentStyle','#');
fclose(file);
phi=C{1};
t=C{2};
ntor=C{3};
npol=C{4};
size=length(t);

% Anzahl Stützstellen von phi
Np=1;
for i=2:size
    if(t(i)==t(i-1)) 
        Np=Np+1;
     else break;
    end
end
Nt=size/Np;    %Anzahl Stützstellen von t

% Matrizen setzen
X=phi(1:Np);
Y=t(1:Np);
Ztor=ntor(1:Np);
Zpol=npol(1:Np);

for i=2:Nt
    X=[X,phi((i-1)*Np+1:i*Np)];
    Y=[Y,t((i-1)*Np+1:i*Np)];
    Ztor=[Ztor,ntor((i-1)*Np+1:i*Np)];
    Zpol=[Zpol,npol((i-1)*Np+1:i*Np)];
end

% Phi Range setzen
imin=fix(0.5*phimin/pi*Np)+1;   % fix rounds towards zero
imax=fix(0.5*phimax/pi*Np);
if(phimin < 0)
    imin=1;
end
X2=X(imin:imax,:);
Y2=Y(imin:imax,:);
Ztor2=Ztor(imin:imax,:);
Zpol2=Zpol(imin:imax,:);
if(phimin < 0)
    phimin=phimin+2*pi;
    imin=fix(0.5*phimin/pi*Np)+1;
    X=[X(imin:Np-1,:)-2*pi;X2];
    Y=[Y(imin:Np-1,:);Y2];
    Ztor=[Ztor(imin:Np-1,:);Ztor2];
    Zpol=[Zpol(imin:Np-1,:);Zpol2];
else
    X=X2;
    Y=Y2;
    Ztor=Ztor2;
    Zpol=Zpol2;
end

% Plot
clf reset;
schrift=18;

% Figure 1: ntor
%figure(1);
clf;
axes('FontName','Times','FontSize',schrift);
contourf(X,Y,Ztor,b);
%h=surf(X,Y,Z);
%shading interp
axis tight;
%ylim([0.195 0.477])
title('n_{tor}','FontName','Times','FontSize',schrift);
xlabel({'j'},'FontName','Symbol','FontSize',schrift);
ylabel({'t'},'FontName','Times','FontSize',schrift);
zlabel('n_{tor}','FontName','Times','FontSize',schrift);
h=findobj('Type','patch'); %nur für contourf
set(h,'EdgeColor','none');
% contourcmap(b,'jet','colorbar','on', ...
%     'FontName','Times','FontSize',schrift);
colormap(jet(length(b)-1));
colorbar('FontName','Times','FontSize',schrift);

F(j)=getframe(gcf);
[G,Map] = frame2im(F(j));
%mov = addframe(mov,F(j));
H(:,:,:,j)=G;
end
mov = immovie(H,Map);
movie2avi(mov,'test.avi')
%movie(F,5);
%mov = close(mov);
