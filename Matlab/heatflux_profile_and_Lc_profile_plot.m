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
yL = yL1;%(yL1+yL2)/2;
yLn = yL;%/max(yL)*2.5;   % con. length in arbitrary units

figure
% plot(x,smooth(y,5),'k-','linewidth',1);
% hold on
% plot(x,yLn,'k--','linewidth',1);
% set(gca,'linewidth',1.5)
% xlim([-1.5 12]);
%grid on
[AX,H1,H2] = plotyy(x,smooth(y,5),x,yLn,'plot');

xlabel('s_{wall} [cm]');

set(get(AX(1),'Ylabel'),'String','q [10^6 W/m^2]','Color','k') 
set(get(AX(2),'Ylabel'),'String','L_c [km]','Color','k','Rotation',-90,'VerticalAlignment','bottom') 
set(H1,'Color','k','LineStyle','-','linewidth',1.5)
set(H2,'Color','k','LineStyle','--','linewidth',1.5)
set(AX(1),'YColor','k','xlim',[-1.5 12],'ylim',[0 2.5],'linewidth',1.5,'box','off','YTick',0:0.5:2.5, 'XTick',[0:2:15])
set(AX(2),'YColor','k','xlim',[-1.5 12],'ylim',[0 1.5],'linewidth',1.5,'box', 'off','XAxisLocation', 'top','XTickLabel',[],'YTick',0:0.25:1.5, 'XTick',[0:2:15])

% ylabel('q [10^6 W/m^2]');
% text(1.05,0.65,'L_c [a.u.]','Rotation',-90, 'Units', 'normalized');

%------------------------------------------------------------------------

load 129194_3000_3200ms_D_alpha_DiMES.txt
x2 = X129194_3000_3200ms_D_alpha_DiMES(:,1)/10.0;
y2 = X129194_3000_3200ms_D_alpha_DiMES(:,3)/100.0;

load 129194_3000_3200ms_heatFlux_IR.txt
x = X129194_3000_3200ms_heatFlux_IR(:,1)/10.0;
y = X129194_3000_3200ms_heatFlux_IR(:,3)/1e+6;

figure
[AX,H1,H2] = plotyy(x,y,x2,y2,'plot');
xlabel('s_{wall} [cm]');
set(get(AX(1),'Ylabel'),'String','q [10^6 W/m^2]','Color','k')
set(get(AX(2),'Ylabel'),'String','D_{\alpha} intensity [a.u.]','Color','k','Rotation',-90,'VerticalAlignment','bottom')
set(H1,'Color','k','LineStyle','-','linewidth',1.5)
set(H2,'Color','k','LineStyle','--','linewidth',1.5)
set(AX(1),'YColor','k','xlim',[-10 20],'ylim',[0 2.5],'linewidth',1.5,'box','off','YTick',0:0.5:2.5, 'XTick',[-8:4:20])
set(AX(2),'YColor','k','xlim',[-10 20],'ylim',[0.5 3],'linewidth',1.5,'box', 'off','XAxisLocation', 'top','XTickLabel',[],'YTick',0.5:0.5:3, 'XTick',[-8:4:20])





