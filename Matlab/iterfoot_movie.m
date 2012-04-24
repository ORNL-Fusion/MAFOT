for i=118:4:118
    Target = 2; % 1=inner 2=outer target
    iterfoot_function(['foot_out_' num2str(i) '-0-' num2str(116) '_rd.dat'],i,Target);
    %close all
end