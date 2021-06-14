function out = get_contour_levels(C,xo,yo)

Npoints = 0;
kk = 0;
ii = 0;
while kk+Npoints < length(C)
    kk = kk + Npoints + 1;
    Npoints = C(2,kk);
    if C(1,kk+1) == C(1,kk+Npoints) && C(2,kk+1) == C(2,kk+Npoints)
        [in,~] = inpolygon(xo,yo,C(1,kk+1:kk+Npoints),C(2,kk+1:kk+Npoints));
        if in == 1
            ii = ii + 1;
            Rb{ii} = C(1,kk+1:kk+Npoints);
            Zb{ii} = C(2,kk++1:kk+Npoints);
            V(ii) = C(1,kk);
            N(ii) = C(2,kk);
        end
    end
end

out.C = C;
out.value = V;
out.n_points = N;
out.Rb = Rb;
out.Zb = Zb;

end
