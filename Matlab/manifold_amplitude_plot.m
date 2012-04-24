for k=56:2:56

load(['manifold_amplitude_' num2str(k) '-0-50_cut']);
Rout_cut = Rout;
Sout_cut = Sout;
load(['manifold_amplitude_' num2str(k) '-0-50_all']);
Rout_all = Rout;
Sout_all = Sout;
load manifold_amplitude_off-0-50_cut

% figure;
% plot(Rout,Sout,'r-',Rout_cut,Sout_cut,'b-',...
%     [0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')

% links und rechts abschneiden
i = min(find(Rout<=5.2));
j = min(find(Rout<=8.15));
Rout = Rout(j:i);
Sout = Sout(j:i);
i = min(find(Rout_cut<=5.2));
j = min(find(Rout_cut<=8.15));
Rout_cut = Rout_cut(j:i);
Sout_cut = Sout_cut(j:i);
i = min(find(Rout_all<=5.2));
j = min(find(Rout_all<=8.15));
Rout_all = Rout_all(j:i);
Sout_all = Sout_all(j:i);

% figure(3);
% plot(Rout,Sout,'r-',Rout_cut,Sout_cut,'b-',...
%     [0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')

% zweideutigen Abschnitt in off-0-50 entfernen
diff = Rout - [10;Rout(1:end-1)];
i = min(find(diff>0))-1;        % 455
j = min(find(Rout<Rout(i)));    % 474
Rout = [Rout(1:i);Rout(j:end)];
Sout = [Sout(1:i);Sout(j:end)];

% evtl. zweideutige Abschnitte in allen anderen entfernen
diff_cut = Rout_cut - [10;Rout_cut(1:end-1)];
while(length(find(diff_cut>0)) > 0)
    i = min(find(diff_cut>0))-1;
    j = min(find(Rout_cut<Rout_cut(i)));
    Rout_cut = [Rout_cut(1:i);Rout_cut(j:end)];
    Sout_cut = [Sout_cut(1:i);Sout_cut(j:end)];
    diff_cut = Rout_cut - [10;Rout_cut(1:end-1)];
end

% figure(4);
% plot(Rout,Sout,'r-',Rout_cut,Sout_cut,'b-',...
%     [0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')

% Interpolieren und addieren
pp = spline(Rout,Sout);
Sout_spline = ppval(pp,Rout_cut);

Sout_sum = Sout_cut + Sout_spline;

% Plot
figure;
plot(Rout_all,Sout_all,'k-',Rout_cut,Sout_sum,'r-',...
    [0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')

xlim([5 8.5]);
ylim([-1.5 1.5]);

if(k<10)
    text(7,-1.2,['Phase: 0' num2str(k) ' deg']);
else
    text(7,-1.2,['Phase: ' num2str(k) ' deg']);
end


filenameout = ['manifold_amplitude_all_' num2str(k/2 + 1)];
print('-r300', '-djpeg', filenameout)

close all
end
