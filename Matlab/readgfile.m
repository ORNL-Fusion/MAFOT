% readgfile.m - global version
%
% First Created: 			2006			Ilon Joseph
% Last Modified: 			2007/08/03	Oliver Schmitz
%function [Rg, Zg, PsiN] = readgfile(shotnum, tM)

%clear
shotnum=132741;
tM=3000;
gfile=['~/c++/d3d/gfiles/g' int2str(shotnum) '.0' int2str(tM)];

global Br
global Bz
global Rg
global Zg
global Rt0
global Bt0
global Raxis
global Zaxis


shot=shotnum;
time=tM;

fid = fopen(gfile,'r');

tag = fscanf(fid,'%s',4);

b= fscanf(fid,'%i',3);
nr=b(2);
nz=b(3);

g = fscanf(fid,'%e',[5,4]);
g=g';
Xdim=g(1,1);
Zdim=g(1,2);
Rt0=g(1,3);      %for bt0
R1=g(1,4);       %inner edge
Zmid=g(1,5);
Raxis=g(2,1);
Zaxis=g(2,2);
PsiAxis=g(2,3);
PsiSep=g(2,4);
Bt0=g(2,5);      %Btor in RHS system
Ip =g(3,1);      %Itor in RHS system

Fpol=fscanf(fid,'%e',nr);
Pres=fscanf(fid,'%e',nr);
FFprime=fscanf(fid,'%e',nr);
Pprime=fscanf(fid,'%e',nr);
PsiRZ=fscanf(fid,'%e',[nr,nz]);
PsiRZ=PsiRZ';
qpsi=fscanf(fid,'%e',nr);

b=fscanf(fid,'%i',2);
Nlcfs=b(1);
Nwall=b(2);
lcfs=fscanf(fid,'%e',2*Nlcfs);
wall=fscanf(fid,'%e',2*Nwall);

fclose(fid);

nz2=floor(nz/2);
[Rg,Zg]=meshgrid( R1 + Xdim*(0:nr-1)/(nr-1), Zmid + Zdim*(-nz2:nz2)/(nz-1)  );
PsiN = (PsiRZ-PsiAxis)/(PsiSep-PsiAxis);
PsiNMax=max(max(PsiN));

figure;
clf;
% %surf(Rg,Zg,-PsiN)
% shading interp
% lighting phong 
% camlight
% view(2)
% 
hold on
contour(Rg,Zg,PsiN,(0:0.1:0.9),'k','Linewidth',1)
contour(Rg,Zg,PsiN,(0.999999:0.000001:1),'r','Linewidth',2)
contour(Rg,Zg,PsiN,1.02:0.02:1.04,'k','Linewidth',1)
% C = contour(Rg,Zg,PsiN,0.95,'k','Linewidth',1);
% 
% % I coils
% rI=[2.184,2.394,2.394,2.184];
% zI=[1.012,0.504,-0.504,-1.012];
% %plot(rI(:),zI(:),'ok','MarkerFaceColor','r','MarkerSize',10);
% rImid=(rI(1)+rI(2))/2;
% zImid=(zI(1)+zI(2))/2;
% 
% % C coils
% rC=[3.23, 3.23];
% zC=[0.8,-0.8];
% %plot(rC(:),zC(:),'ok','MarkerFaceColor','b','MarkerSize',10)
% 
rwall=wall(1:2:2*Nwall)';
zwall=wall(2:2:2*Nwall)';
rlcfs = lcfs(1:2:2*Nlcfs)';
zlcfs = lcfs(2:2:2*Nlcfs)';
plot(rlcfs,zlcfs,'r','Linewidth',2)
plot(rwall,zwall,'k','Linewidth',2)
plot(Raxis,Zaxis,'k+')
axis equal tight
% 
% title(horzcat(num2str(shot),'.0',num2str(time)))
% set(get(gca,'title'),'fontsize',14,'fontweight', 'bold')
% 
% 
% %Br = -cdiff2(Zg,PsiRZ);
% %Bz = cdiff2(Rg',PsiRZ');
% %Bz = Bz';
% %Bt = Bt0*Rt0./Rg;
% %figure(2)
% %s=4
% %quiver(Rg(1:s:129,1:s:129),Zg(1:s:129,1:s:129),...
% %    Br(1:s:129,1:s:129),Bz(1:s:129,1:s:129),6)
% %axis equal tight
