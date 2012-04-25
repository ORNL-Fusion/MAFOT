function [R Z] = t2toroidal(t,target)
%[R Z] = t2toroidal(t,target)
%target = 0 (inner target)
%target = 1 (outer target)

    % Wall parameter
    WallFile = '/home/rack/JET/wall.dat';
    if (target==0)
        target_start = 97;
        target_end = 80;
    else
        target_start = 97;
        target_end = 159;
    end
    wallfile=fopen(WallFile);
    C=textscan(wallfile,'%f%f','commentStyle','#');
    fclose(wallfile);
    Rwall=C{1};
    Zwall=C{2};
    target_max = length(Rwall);

    if (target_start > target_end)
        Rdiv = Rwall(target_start:-1:target_end);
        Zdiv = Zwall(target_start:-1:target_end);
    else
        if (target_end > target_max)
            Rdiv = [Rwall(target_start:target_max); Rwall(2:(target_end-target_max+1))];
            Zdiv = [Zwall(target_start:target_max); Zwall(2:(target_end-target_max+1))];
        else
            Rdiv = Rwall(target_start:target_end);
            Zdiv = Zwall(target_start:target_end);
        end
    end
    
    N = length(Rdiv)-1;
    l = zeros(N,1);
    
    for k=1:N
        l(k) = sqrt((Rdiv(k+1)-Rdiv(k))^2 + (Zdiv(k+1)-Zdiv(k))^2)*100;
    end

    %disp(['Der innere Divertor hat ' num2str(sum(l)) ' cm']);
    %disp(['Der äußere Divertor hat ' num2str(sum(l)) ' cm']);

    if (t<0), disp('t must be positive'); end
    %-------Finden des entsprechenden Gradenstücks-------------
    k = 1;
    while (t > 0)
        if (k > N)
            disp(['t is higher then the allowed value of ' num2str(sum(l)) ' cm'])
            return;
        end
        t = t-l(k);
        k = k+1;
    end
    %------Berechnung der Geradengleichung dieses Gradenstücks-
    if (t==0)
        R = Rdiv(k);
        Z = Zdiv(k);
    else
        t = t/l(k-1);
        R = (Rdiv(k)-Rdiv(k-1))*t + Rdiv(k);
        Z = (Zdiv(k)-Zdiv(k-1))*t + Zdiv(k);
    end
end