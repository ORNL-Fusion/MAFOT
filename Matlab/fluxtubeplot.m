clear;

% Anzahl der Streamlines in File
N = 4*7;  

% Filename
file=fopen('struct_up.dat');
C=textscan(file,'%f%f%f%f%f','commentStyle','#');
fclose(file);
x = C{1};
y = C{2};
z = C{3};
r = C{4};
phi = C{5};

file=fopen('struct_right.dat');
C=textscan(file,'%f%f%f%f%f','commentStyle','#');
fclose(file);
x = [x;C{1}];
y = [y;C{2}];
z = [z;C{3}];
r = [r;C{4}];
phi = [phi;C{5}];

file=fopen('struct_down.dat');
C=textscan(file,'%f%f%f%f%f','commentStyle','#');
fclose(file);
x = [x;C{1}];
y = [y;C{2}];
z = [z;C{3}];
r = [r;C{4}];
phi = [phi;C{5}];

file=fopen('struct_left.dat');
C=textscan(file,'%f%f%f%f%f','commentStyle','#');
fclose(file);
x = [x;C{1}];
y = [y;C{2}];
z = [z;C{3}];
r = [r;C{4}];
phi = [phi;C{5}];

begin = zeros(1,N)';
ende = zeros(1,N)';

j = 2;
begin(1) = 1;
for i = 1:length(phi)-1
    if(phi(i)>0 && phi(i+1)<0) 
        begin(j) = i+1;
        j = j+1;
    end
end
for i = 1:N-1
    ende(i) = begin(i+1)-1;
end
ende(N) = length(phi);
center = find(phi==0);
extend = ende - begin + 1;

maxL = round(max([18*phi(ende)/pi;-18*phi(begin)/pi]));
minL = round(min([18*phi(ende)/pi;-18*phi(begin)/pi]));

R = NaN*ones(2*maxL+1,N);
PHI = NaN*ones(2*maxL+1,N);
Z = NaN*ones(2*maxL+1,N);
X = NaN*ones(2*maxL+1,N);
Y = NaN*ones(2*maxL+1,N);

shift = maxL+1;
for j=1:N
    R(shift,j) = r(center(j));
    PHI(shift,j) = phi(center(j));
    Z(shift,j) = z(center(j));
    
    X(shift,j) = x(center(j));
    Y(shift,j) = y(center(j));
end

for i=1:minL
    for j=1:N
        if(center(j)+i<=ende(j))    R(shift+i,j)   = r(center(j)+i);   end
        if(center(j)+i<=ende(j))    PHI(shift+i,j) = phi(center(j)+i); end
        if(center(j)+i<=ende(j))    Z(shift+i,j)   = z(center(j)+i);   end
        if(center(j)-i>=begin(j))   R(shift-i,j)   = r(center(j)-i);   end
        if(center(j)-i>=begin(j))   PHI(shift-i,j) = phi(center(j)-i); end
        if(center(j)-i>=begin(j))   Z(shift-i,j)   = z(center(j)-i);   end
        
        if(center(j)+i<=ende(j))    X(shift+i,j)   = x(center(j)+i);   end
        if(center(j)+i<=ende(j))    Y(shift+i,j)   = y(center(j)+i);   end
        if(center(j)-i>=begin(j))   X(shift-i,j)   = x(center(j)-i);   end
        if(center(j)-i>=begin(j))   Y(shift-i,j)   = y(center(j)-i);   end

    end
end

%h = surf(PHI,R,Z);
%h = surf(X,Y,Z);
%set(h,'EdgeColor','none');
%mesh(X,Y,Z);
%mesh(PHI,R,Z);

figure(1);
h = surfl(X,Y,Z);
set(h,'EdgeColor',[.5 .5 .5])
colormap(hsv);
axis vis3d
view(-37,30);
%axis off
light('Position',[2 -4 5])
light
%hold off

figure(2);
wert = 4.5*36*pi/18
wo = find(round(1000*(PHI-wert))/1000==0)
plot(R(wo),Z(wo),'kx-','markerSize',3);