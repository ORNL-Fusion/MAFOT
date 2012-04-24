
FileToOpen = 'plot_all_large.dat';
itt = 300;
N = 120;

file=fopen(FileToOpen);
C=textscan(file,'%f%f%f%f%f%f','commentStyle','#');
fclose(file);
theta=C{1};
r=C{2};
phi=C{3};
psi=C{4};
R=C{5};
Z=C{6};


for i = 1:N
    index = itt*i;
    x = phi(index);
    if(x < itt*360)
        psi((i-1)*itt+1)
        break;
    end    
end