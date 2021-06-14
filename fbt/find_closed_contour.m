function out = find_closed_contour(C,params)
%% C is the output of the function contourc

out.rb = [];
out.zb = [];
k = 1;
n = size(C,2);
while k < n
    nk = C(2,k);
    if C(1,k+1) == C(1,k+nk) && C(2,k+1) == C(2,k+nk) && isempty(find(inpolygon(C(1,k+1:k+nk),C(2,k+1:k+nk),params.limiter(:,1),params.limiter(:,2)) == 0,1))
        % Closed contour
        out.rb = C(1,k+1:k+nk);
        out.zb = C(2,k+1:k+nk);
        break
    end
    k=k+nk+1;
end
    
end