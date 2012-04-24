
for i=0:2:88

load(['manifold_amplitude_' num2str(i) '-0-50_cut']);
Rout_cut = Rout;
Sout_cut = Sout;
load manifold_amplitude_off-0-50_cut

figure;
plot(Rout,Sout,'r-',Rout_cut,Sout_cut,'b-',...
    [0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')

xlim([5 8.5]);
ylim([-1.5 1.5]);

if(i<10)
    text(7,-1.2,['Phase: 0' num2str(i) ' deg']);
else
    text(7,-1.2,['Phase: ' num2str(i) ' deg']);
end


filenameout = ['manifold_amplitude_' num2str(i/2 + 1)];
print('-r300', '-djpeg', filenameout)

close all
end