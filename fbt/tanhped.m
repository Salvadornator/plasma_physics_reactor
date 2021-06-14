function y = tanhped(x,ysep,xmid,d,a0,a1,alfa1,alfa2)

xped = xmid - d/2;
y = ysep + a0/2*(1 - tanh(2*(x - xmid)/d)) + a1*heaviside(1 - x/xped).*(1 - (x/xped).^alfa1).^alfa2;

end