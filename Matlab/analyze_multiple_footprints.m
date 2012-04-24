% set printme = 0 in iterfoot_function
clear;
clc;

n = 29;
x = [2:4:114]';

Area = zeros(n,1);
Tip = zeros(n,1);
psi_min = zeros(n,1);
flf98 = zeros(n,1);
flf95 = zeros(n,1);
Lc_average = zeros(n,1);
flf440 = zeros(n,1);
flf1 = zeros(n,1);

for i = 1:n
    FileToOpen = ['foot_out_' num2str(4*i-2) '-0-' num2str(114-4*i+2) '_rd.dat'];
    Target = 2; % 1=inner 2=outer target
    
    [Z,X,Y,Nt,Np] = iterfoot_function(FileToOpen,i,Target);   
    [Area(i),Tip(i),Sp] = iterfootprint_area_function(X,Y,Z,Nt,Np,Target);
    [psi_min(i),flf98(i),flf95(i),Lc_average(i),flf440(i),flf1(i)] = iterfieldlinefraction_function(FileToOpen);
    close all
    disp(i);
end
%Areanorm  = Area/max(Area);
output = [x Area Tip psi_min flf98 flf95 Lc_average flf440 flf1];
if(Target == 2), save('all_data_out_rd.dat', 'output','-ascii');
else save('all_data_in.dat', 'output','-ascii');
end