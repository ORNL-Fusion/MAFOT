function [X,Y,Z] = read_foot_lam_data(FileToOpen,WhatShallIPlot,yVariedFirst)
%[X,Y,Z] = read_foot_lam_data(FileToOpen,WhatShallIPlot,yVariedFirst)
%reads footprint and laminar data and returns X,Y,Z 
%X,Y,Z are given as stored in file and maybe must be modified for
%differenz plots

    if (exist('yVariedFirst','var')==0), yVariedFirst=0; end
    % --- Read data from input file -----------------------------------
    file=fopen(FileToOpen);
    % die Matrizen hier sind transponiert zu den unteren, ist aber egal
    %if(WhatShallIPlot==2)
        C=textscan(file,'%f%f%f%f%f','commentStyle','#');
    %else
    %    C=textscan(file,'%f%f%f%f%*[^\n]','commentStyle','#');
    %end
    fclose(file);
    phi=C{1};
    t=C{2};
    val=C{WhatShallIPlot+3};
    sizeofvec=length(t);

    % # grid points of phi
    Np=1;
    for i=2:sizeofvec
        if(t(i)==t(i-1)) 
            Np=Np+1;
        else break;
        end
    end
    Nt=sizeofvec/Np;    % # grid points of t

    if(yVariedFirst==1)
        Nt=1;
        for i=2:sizeofvec
            if(phi(i)==phi(i-1)) 
                Nt=Nt+1;
            else break;
            end
        end
        Np=sizeofvec/Nt;    % reset # grid points of phi
    end

    X = reshape(phi,Np,Nt);
    Y = reshape(t,Np,Nt);
    Z = reshape(val,Np,Nt);
end