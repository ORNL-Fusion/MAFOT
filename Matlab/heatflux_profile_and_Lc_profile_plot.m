load con_length_150deg_129194_3000_inner_fl.dat
load con_length_150deg_129194_3000_inner_ip.dat
load heatflux_150deg_129194_3000_inner_ip.dat

x = heatflux_150deg_129194_3000_inner_ip(:,1);
y1 = heatflux_150deg_129194_3000_inner_ip(:,2);
% y1 = Mittlere Energie pro Teilchen in [keV]
% W‰rmefluﬂ = 1/2 Teilchendichte * thermische Geschwindigkeit * Energie pro Teilchen
% q = 1/2 n_i * sqrt(kT/m_i) * E        mit [E] = J
%   = 7.83e-16 * n_i * T^(3/2) W/m^2    mit [T] = eV
% Annahme: n_i = 3.6e19 1/m^3
y = 28188*(y1*1000).^1.5/1e6;   % in 1e6 W/m^2

yL1 = con_length_150deg_129194_3000_inner_fl(:,2);
yL2 = con_length_150deg_129194_3000_inner_ip(:,2);
yL = (yL1+yL2)/2;
yLn = yL/max(yL)*2.5;   % con. length in arbitrary units

figure
plot(x,yLn,'b-',x,y,'r-')
xlim([-5 15]);
grid on

xlabel('s_{wall} [cm]');
ylabel('q [10^6 W/m^2] (red line)');
text(1.05,0.65,'L_c [a.u.] (blue line)','Rotation',-90, 'Units', 'normalized');

