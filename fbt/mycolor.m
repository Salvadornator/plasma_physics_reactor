function c = mycolor

c = zeros(12,3);
c(1,:) = [237 32 36];   % Red
c(2,:) = [164 194 219]; % Light blue
c(3,:) = [57 83 164];   % Dark Blue
c(4,:) = [50 180 80];   % Light Green
c(5,:) = [0 0 0];       % Black
c(6,:) = [0 110 0];     % Green
c(7,:) = [256 0 256];   % Magenta
c(8,:) = [0 256 256];   % Cyan
c(9,:) = [0 0 256];     % Normal blue
c(10,:) = [256 0 0];    % Normal red 
c(12,:) = [256 102 154];% Lilac
c = c/256;

end