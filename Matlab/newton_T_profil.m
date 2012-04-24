% Newton für heatflux_footprint_generate.m

function x = newton_T_profil(x0,wert)
epsilon = 1e-14;
dx = inf;
n = 0;

while (max(max(abs(dx))) > epsilon)
    [y, dy] = T_profil(x0);
    dx = (y-wert)./dy;
    x0 = x0 - dx;
    n = n + 1;
    if(n > 20) disp('Newton: Keine Konvergenz'); x0 = NaN; break; end
end

x = x0;
end

%--- End Main -------------------------------------------------------------
%--------------------------------------------------------------------------

%--- T_profil -------------------------------------------------------------
function [T, dT] = T_profil(psi)
xs = 0.975;     % Symmetry point in Pedestal
dw = 0.04;      % half of Pedestal width
Tmax = 2;       % Temparature within main plasma

T = 0.5*tanh(2*(xs-psi)/dw) + 0.5;  % normiert auf 1 im Inneren
dT = -1/dw*(1-(2*T-1)*(2*T-1));
T = Tmax*T;    dT = Tmax*dT;              % T im Inneren jetzt 2 keV
end
