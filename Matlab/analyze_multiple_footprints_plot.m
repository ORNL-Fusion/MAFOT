clear
close all

% load all_data_out_ok.dat
% data = all_data_out_ok;
% 
% Area1 = data(:,2);
% Tip1 = data(:,3);
% psi_min1 = data(:,4);
% flf981 = data(:,5);
% flf951 = data(:,6);
% Lc_average1 = data(:,7);
% flf4401 = data(:,8);
% flf11 = data(:,9);

load all_data_out.dat
data = all_data_out;

x = data(:,1);
Area = data(:,2);
Tip = data(:,3);
psi_min = data(:,4);
flf98 = data(:,5);
flf95 = data(:,6);
Lc_average = data(:,7);
flf440 = data(:,8);
flf1 = data(:,9);


% Area
figure
plot(x,Area/max(Area),'k-')%,x,Area1/max(Area1),'k--')
xlabel('upper coil Phase [deg]');
ylabel('Footprint Area (a.u.)');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% Tip
figure
plot(x,Tip,'k-')%,x,Tip1,'k--')
xlabel('upper coil Phase [deg]');
ylabel('Distance: Tip <-> SL [cm]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% psi_min
figure
plot(x,psi_min,'k-')%,x,psi_min1,'k--')
xlabel('upper coil Phase [deg]');
ylabel('\fontname{Symbol}y _{\fontname{Times}Min}');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% flf98
figure
plot(x,flf98,'k-')%,x,flf981,'k--')
xlabel('upper coil Phase [deg]');
ylabel('flf_{98} [%]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% flf95
figure
plot(x,flf95,'k-')%,x,flf951,'k--')
xlabel('upper coil Phase [deg]');
ylabel('flf_{95} [%]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% Lc
figure
plot(x,Lc_average,'k-')%,x,Lc_average1,'k--')
xlabel('upper coil Phase [deg]');
ylabel('average L_c [km]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% flf440
figure
plot(x,flf440,'k-')%,x,flf4401,'k--')
xlabel('upper coil Phase [deg]');
ylabel('flf_{440m} [%]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on

% flf1k
figure
plot(x,flf1,'k-')%,x,flf11,'k--')
xlabel('upper coil Phase [deg]');
ylabel('flf_{1km} [%]');
set(gca,'XTick', 0:10:10*ceil(max(x/10)));
grid on