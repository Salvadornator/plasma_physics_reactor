function k = imin(varargin)

% IMIN	Index of minimum
%    IMIN returns the index of minimum (second output of MIN)

[tmp,k] = min(varargin{:});
