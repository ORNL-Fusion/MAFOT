clear;

% Parameter
jmax=36;    % Anzahl der Files für Video
direction=-1;   % phi Richtung: 1=pos   -1=neg

% Namen
nameni='d3dman_2690fix_stl1_I_';
namenc='d3dman_2690fix_stl1_C_';
namenic='d3dman_2690fix_stl1_I+C_';

scrsz = get(0,'ScreenSize');
figure('Position',[0.1*scrsz(3),0.15*scrsz(4),0.8*scrsz(3),0.7*scrsz(4)])

% Ungestörte Separatrix
file=fopen('d3dman_2690fix_stl1_noRMP_0_RZ.dat');
C=textscan(file,'%f%f%f','commentStyle','#');
fclose(file);
R=C{1};
Z=C{2};

for j=[1 1:jmax]    % 1. Frame wird zweimal verarbeitet, das stabilisiert den ersten Frame
    % Einlesen der Daten aus Main File
    if(direction==-1) num=int2str((jmax-j+1)*10);
    else num=int2str(j*10);
    end
    name=[nameni num '_0_RZ.dat'];
    filei=fopen(name);
    C=textscan(filei,'%f%f%f','commentStyle','#');
    fclose(filei);
    Ri=C{1};
    Zi=C{2};

    name=[namenc num '_0_RZ.dat'];
    filec=fopen(name);
    C=textscan(filec,'%f%f%f','commentStyle','#');
    fclose(filec);
    Rc=C{1};
    Zc=C{2};

    name=[namenic num '_0_RZ.dat'];
    fileic=fopen(name);
    C=textscan(fileic,'%f%f%f','commentStyle','#');
    fclose(fileic);
    Ric=C{1};
    Zic=C{2};

    % Plot
    clf reset;
    %clf;
    schrift=18;
    axes('FontName','Times','FontSize',schrift);
    hold on
    plot(R,Z,'k-','LineWidth',2);
%     plot(Ri,Zi,'r-','LineWidth',2);
%     plot(Rc,Zc,'g-','LineWidth',2);
    plot(Ric,Zic,'b-','LineWidth',2);
    hold off
    legend('noRMP','I','C','I+C');
    axis tight;
    xlim([1.179 1.429]);
    ylim([-1.47 -0.85]);
    xlabel({'R [m]'},'FontName','Times','FontSize',schrift);
    ylabel({'Z [m]'},'FontName','Times','FontSize',schrift);

    % Get Movie Frame
    F(j)=getframe(gcf);
    [G,Map] = frame2im(F(j));
    %mov = addframe(mov,F(j));
    H(:,:,:,j)=G;
end

% Build Movie
mov = immovie(H,Map);
movie2avi(mov,'2690_manifolds_I+C.avi','quality',100)
%movie(F,5);
%mov = close(mov);
