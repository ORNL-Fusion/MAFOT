for j=56:2:56
    FileToOpen = ['man_stl1_' num2str(j) '-0-50_all_0.dat'];
%FileToOpen = 'man_stl1_56-0-50_all_0.dat';

    file = fopen(FileToOpen);
    C = textscan(file,'%f%f%f%f%f','commentStyle','#');
    fclose(file);
    R = C{1};
    Z = C{2};

%     figure
%     plot(R,Z,'k-')

    N = length(R);
    Rout = zeros(N,1);
    Sout = zeros(N,1);
    load C:\C++\ITER\800001.00000\Sep;
    skipit=1;
    pp = spline(Sep(1,:),Sep(2,:));

    for i=2:N
        if((R(i)<8.19 | Z(i)>0.4) & skipit==1) continue; 
        else skipit=0;
        end

        fun = @(x) (x-R(i))^2+(ppval(pp,x)-Z(i))^2;
        x0 = Rout(i-1);
        options = optimset('TolFun',1e-4);
        [x,fval] = fminsearch(fun,x0,options);
        Rout(i) = x;
        y = [x,ppval(pp,x)]-[R(i),Z(i)];
        Sout(i) = sign(y(2))*norm(y);
        if(mod(N-i,100)==0) 
            clc
            disp([1 j N-i]); 
        end
    end
    Sout=Sout(find(Rout~=0));
    Rout=Rout(find(Rout~=0));
    
    figure
    plot(Rout,Sout)
    
    filetosave = ['manifold_amplitude_' num2str(j) '-0-50_all'];
%    filetosave = 'manifold_amplitude_off-0-50_all';
    save(filetosave,'Rout','Sout');
    %plot(Rout,Sout,'k-',[0.95*min(Rout) 1.05*max(Rout)],[0 0],'k--')
end

