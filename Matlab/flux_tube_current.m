clear
%------------ Set parameter here --------------------------
% Length [m]
L = 55.3;

% norm. Flux psi 
psi = 0.9987;  

% Area [cm^2]
A = 14.05;   

%----------------------------------------------------------
%----------------------------------------------------------
% Code calculates Current in flux tube 
% Scaling Law: I = j0*A*L0/L*[T(psi)/T(psi0)]^(3/2)
% Parameters:
%               j0 = 125 A/cm^2
%               L0 = 50 m
%               psi0 = 1
%               T(psi0) = 125.5 eV  <-  based on Te profile
% 
% Uses Te profile from 132741@3000ms, 
% fittet with exponential decay between psi = 0.96985 -> 1.04544
% Te = exp(69.91538-72.21462*psi) + 0.02534

% Code Returns:
%               Te [keV]
%               I [A]
%----------------------------------------------------------
%----------------------------------------------------------
% Scaling parameter
j0 = 125;
L0 = 50;
psi0 = 1;

% Te profile
Te0 = exp(69.91538-72.21462*psi0) + 0.02534;
Te = exp(69.91538-72.21462*psi) + 0.02534

% Current
I = j0*A*L0/L*(Te/Te0)^(3/2)
