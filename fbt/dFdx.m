function f = dFdx(x,F,varargin)

if length(x) == length(F)
    Lx = length(x);
    f = zeros(size(x));
    if nargin > 2
        iorder = varargin{1};
    else
        iorder = 1;
    end
    if iorder == 1
        for ii = 2:Lx-1
            f(ii) = (F(ii+1) - F(ii-1))/(x(ii+1) - x(ii-1));
        end
        f = interp1(x(2:Lx-1),f(2:Lx-1),x,'pchip');
    elseif iorder == 2
        for ii = 3:Lx-2
            f(ii) = (-F(ii + 2) + 8*F(ii+1) - 8*F(ii-1) + F(ii-2))/6/(x(ii+1) - x(ii-1));
        end
        f = interp1(x(3:Lx-2),f(3:Lx-2),x,'pchip');
    end
else
    error('Vectors must have same dimension')
end

end