function [RZ] = s2toroidal(s,RorZ)
%[RZ] = s2RZ(s,RorZ)
%convert s to R or Z depending on parameter RorZ
%RorZ, 1: R; 2: Z
%s have to be increasing

    % read converter data
    file = fopen('/home/rack/JET/srz_2006.dat');
    C = textscan(file,'%f%f%f','commentStyle','#');
    fclose(file);
    conv_data = [C{1} C{2} C{3}];
    
    N = length(s);
    RZ = zeros(N,1);
    
    for k = 1:N
        idx = find(s(k) >= conv_data(:,1), 1, 'last');
        if (s(k) == conv_data(idx,1))
            RZ(k) = conv_data(idx,RorZ+1);
        else
            RZ(k) = conv_data(idx,RorZ+1) + ...
                (conv_data(idx+1,RorZ+1) - conv_data(idx,RorZ+1)) / ...
                (conv_data(idx+1,1) - conv_data(idx,1)) * ...
                (s(k) - conv_data(idx,1));
        end
    end
end