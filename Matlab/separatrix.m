clear
%------------ Set parameter here --------------------------
% Shot
shot = 132741;
time = 2630;

% Point position
R = 1.03657259783368;  Z = -1.24943386148447;   

% Hide figures of EFIT and separatrix
close_figures = 0;   % 0: no   1: yes

%----------------------------------------------------------
%----------------------------------------------------------
% Code calculates normalized flux and 
% distance from Separatrix for any Point (R,Z) 
% in the poloidal cross-section
% 
% Inner and outer Strike-Point positions are calculated
% also Distance from inner SP to Target-Edge 
%           -> transform t to Distance from SP
%
% Uses g-file
%----------------------------------------------------------
%----------------------------------------------------------

% Read g-file
[Rg, Zg, PsiN] = readgfile(shot,time);

% Interpolate PsiN(R,Z)
psi = interp2(Rg,Zg,PsiN,R,Z)

% All contour lines with psi = 1
C = contourc(Rg(1,:),Zg(:,1),PsiN,[1 1]);

% Get Separatrix from all contour lines in C <- longest contour line
index = find(C(2,:)==round(C(2,:)) & C(2,:)>3);
max_number = max(C(2,index));
max_index = index(C(2,index)==max_number);
sep = C(:,max_index+1:1:max_index+max_number);

% Add second longest contour line
if(sep(1,1) > 1.0161 & sep(2,1) > -1.3664)  % first point in sep is already inside the vessel boundary
    disp('I did it')
    max_number2 = max(C(2,find(C(2,index)~=max_number)));
    max_index2 = index(C(2,index)==max_number2);
    sep = [C(:,max_index2+1:1:max_index2+max_number2), sep];
end

% Plot Separatrix: Check if code works correct
figure;
hold on;
plot(R,Z,'r+','MarkerSize', 8);
plot(sep(1,:),sep(2,:),'k-','LineWidth',2);

% Plot Wall
wallfile=fopen('C:\c++\d3d\wall.dat');
C3=textscan(wallfile,'%f%f%f%f','commentStyle','#');
fclose(wallfile);
phiout=C3{3};
tout=C3{4};
plot(phiout,tout,'k--','LineWidth',2);
axis equal tight

% Calculate minimum distance (in [m]) from Point to Separatrix
min_dist = 100;
for i=1:size(sep,2)
    distance = norm(sep(:,i)-[R;Z]);
    if(distance < min_dist) min_dist = distance; end
end

min_dist

if(close_figures == 1) close all; end

% Get Strike-Point position
fertig = 0;
for i = 2:size(sep,2)
    chk = 0;
    R = sep(1,i);
    Z = sep(2,i);
    if(R > 1.0161 & Z >= -1.22884 & fertig == 0) chk = 1; end	% separatrix enters through center post wall
    if(R >= 1.15285 & Z <= -1.3664 & chk == 0 & fertig < 2) chk = 2; end	% separatrix leaves through outer target

    % 45° Target-Plate
    R1=1.0161;	Z1=-1.22884;    % upper left Point
    R2=1.15285;	Z2=-1.3664;     % lower right Point

    a = (Z2-Z1)/(R2-R1);
    b = Z1 - a*R1;
    Z_R = a*R + b;
    if(Z > Z_R & chk == 0 & fertig == 0) chk = 3; end	% separatrix enters through 45°target
    
    % After wall intersection: interpolate point of intersection
    if(chk > 0)
        R_alt = sep(1,i-1);
        Z_alt = sep(2,i-1);
        dR = R - R_alt;
        dZ = Z - Z_alt;

        if(chk == 1)	% center post wall
            R = 1.0161;
            x = (R - R_alt)/dR;
            Z = Z_alt + x*dZ;
            disp('Inner SP:     R =                          Z =')
            disp([R Z])
            SPdist = -norm([R-R1, Z-Z1])*100
            fertig = 1;
        elseif(chk == 2)	% outer target
            Z = -1.3664;
            x = (Z - Z_alt)/dZ;
            R = R_alt + x*dR;
            disp('Outer SP:     R =                          Z =')
            disp([R Z])
            fertig = 2;
        elseif(chk == 3)	% inner target (45°)
            x = (Z_alt - a*R_alt - b)/(a*dR - dZ);
            R = R_alt + x*dR;
            Z = Z_alt + x*dZ;
            disp('Inner SP:     R =                          Z =')
            disp([R Z])
            SPdist = norm([R-R1, Z-Z1])*100
            fertig = 1;
        end
    end
end

