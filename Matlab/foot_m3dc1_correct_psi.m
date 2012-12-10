val = C{4};   % Lc
L = reshape(val,Np,Nt);
if(Laminar == 0 && Machine == 1)
    L = L(end:-1:1,:);
end


% Minimal con. length to consider
cmin = 0.075;

Z(L <= cmin) = 1;




% if(Physical == 1)
%     SP = min(Y(L > cmin))   % Z in m
%     Tip = max(Y(L > cmin))  % Z in m
%     Extent = abs(SP - Tip)  % in cm
% else
%     SP = max(Y(L > cmin))   % t in cm; for Z in m swap SP and Tip 
%     Tip = min(Y(L > cmin))  % t in cm; for Z in m swap SP and Tip 
%     Extent = abs(SP - Tip)  % in cm
% end
% 
% Z(Z > 0.999) = 1;
% Z(Y < SP) = 1;
% Z(Y > Tip) = 1;




% inside = 0;
% for j=1:Np
%     for i=1:Nt
%         if(Z(j,i) == 1 && inside == 1), inside = 0; end
%         if(Z(j,i) == 1), continue; end
%         if(Z(j,i) < 1 && inside == 1), continue; end
%         if(Z(j,i) < 1 && Z(j,i-1) == 1 && Z(j,i) > 0.988 && inside == 0)
%             %display([j, i])
%             Z(j,i) = 1; 
%             continue;
%         end
%         if(Z(j,i) <= 0.988 && inside == 0), inside = 1; end
%     end
% end
