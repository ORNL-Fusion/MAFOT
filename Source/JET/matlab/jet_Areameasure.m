clear;
old_path = path;
path(old_path,'./functions');
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Code calculates size of fluxtube area
%
% the red area is the measured one
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --- Set parameter here --------------------------------------------------
% Main Filename
FileToOpen = '/.netmount/storage_blast/rack/079459.55918/foot_out_ef.dat';

% Starting point in non-Machine Mode (3 digits behind the dot!!!):
x =  5.112;
y = 18.690;

% Interval of (short) connection length
ValMin = 0.16;
ValMax = 0.272;

% extended options
Periodic = 0;      % 0: plot treat as not periodic    1: plot is periodic
Laminar = -1;      % 0: Footprint  1: Laminar Plot   -1: auto mode
Target = -1;       % 0: Inner  1: Outer  -1: auto mode

% automatic mode (read filename and try to identify the data type)
[Laminar, Target] = ident_data_type(FileToOpen,Laminar, Target);

% Load file
[X,Y,Z] = read_foot_lam_data(FileToOpen,1);

% Search and mark
Z((Z>ValMax) | (Z<ValMin)) = 0;
Z(Z~=0) = 1;

% Preparation for starting points/values
[idr idc] = find((round(X*1000)/1000==x) & (round(Y*1000)/1000==y));
H = [idr idc];

% Algorithem for area detection
Nc = length(Z(1,:));
Nr = length(Z(:,1));
flaeche = 2;
while (isempty(H) == 0)
    idr = H(end,1);
    idc = H(end,2);
    
    H = H(1:end-1,1:2);
    Z(idr,idc) = flaeche;
    
    idcp = idc+1;
    if ((idcp > Nc) && Periodic), idcp = 1; end
    if (idcp <= Nc)
        if (Z(idr,idcp) == 1)
            H(end+1,1:2) = [idr,idcp];
            Z(idr,idcp) = -1;
        end
    end
    idrp = idr+1;
    if ((idrp > Nr) && Periodic), idrp = 1; end
    if (idrp <= Nr)
        if (Z(idrp,idc) == 1)
            H(end+1,1:2) = [idrp,idc];
            Z(idrp,idc) = -1;
        end
    end
    idcm = idc-1;
    if ((idcm < 1) && Periodic), idcm = Nc; end
    if (idcm >= 1)
        if (Z(idr,idcm) == 1)
            H(end+1,1:2) = [idr,idcm];
            Z(idr,idcm) = -1;
        end
    end
    idrm = idr-1;
    if ((idrm < 1) && Periodic), idrm = Nr; end
    if (idrm >= 1)
        if (Z(idrm,idc) == 1)
            H(end+1,1:2) = [idrm,idc];
            Z(idrm,idc) = -1;
        end
    end
end

figure;
clf;
pcolor(X,Y,Z);
shading interp;

% Remove background in the plot
MyColorMap=get(gcf,'Colormap');
MyColorMap(1,:) = 1;
set(gcf,'Colormap',MyColorMap)

% Output
IDs = find(Z==flaeche);
dY = abs(Y(1,2)-Y(1,1));
dX = abs(X(2,1)-X(1,1));
if (Laminar == 0)
    R = zeros(length(IDs),1);
    for k=1:length(IDs)
        [R(k) ~] = t2toroidal(Y(IDs(k)),Target);
    end
    area = sum(dY*dX*R*100);
else
    area = dY*100*dX*100*length(IDs);    
end
%2*pi*275*10*length(find(Z==flaeche))/Nr/Nc
disp(['Area has a size of ' num2str(area) ' cm^2'])

path(old_path);