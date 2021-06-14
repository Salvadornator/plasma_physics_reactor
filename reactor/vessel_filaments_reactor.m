function [W,dRv,dZv] = vessel_filaments_reactor

%%
Wall = load('/home/canal/matlab/equil_reconst/reactor/reactor_vacuum_vessel');
ind = find(isnan(Wall(:,1)));
Rv = [(Wall(1,1) + Wall(ind(1)+1,1))/2 (Wall(1,1) + Wall(ind(1)+1,1))/2 (Wall(3,1) + Wall(ind(1)+3,1))/2 (Wall(3,1) + Wall(ind(1)+3,1))/2 (Wall(1,1) + Wall(ind(1)+1,1))/2];
Zv = [(Wall(1,2) + Wall(ind(1)+1,2))/2 (Wall(2,2) + Wall(ind(1)+2,2))/2 (Wall(3,2) + Wall(ind(1)+3,2))/2 (Wall(4,2) + Wall(ind(1)+4,2))/2 (Wall(1,2) + Wall(ind(1)+1,2))/2];
lv = Zv(2)-Zv(1);
Nv = 55;

% Inner wall
Zvvi = linspace(Zv(1),Zv(2),round(Nv))';
Rvvi = Rv(1)*ones(size(Zvvi));
dlv = mean(diff(Zvvi));
dRin = (Wall(ind(1)+1,1)-Wall(1,1))*ones(size(Rvvi));
dZin = mean(diff(Zvvi))*ones(size(Zvvi));

% Outer wall
Zvvo = linspace(Zv(3),Zv(4),round(Nv))';
Rvvo = Rv(3)*ones(size(Zvvo));
dRout = (Wall(3,1)-Wall(ind(1)+3,1))*ones(size(Rvvo));
dZout = mean(diff(flip(Zvvo)))*ones(size(Zvvo));

% Top wall
lh = Rvvo(1)-Rvvi(1)-5e-3;
Nh = lh/dlv;
Rvvt = linspace(Rvvi(1)+2.5e-3,Rvvo(1)-2.5e-3,round(Nh)+1)';
dlh = mean(diff(Rvvt));
Rvvt = Rvvt(1:end-1)+dlh/2;
Zvvt = linspace(Zv(2),Zv(3),length(Rvvt));
dRtop = mean(diff(Rvvt))*ones(size(Rvvt));
dZtop = (Wall(2,2)-Wall(ind(1)+2,2))*ones(size(Zvvt));

% Bottom wall
Rvvb = Rvvt(end:-1:1);
Zvvb = linspace(Zv(4),Zv(5),length(Rvvt));
dRbot = mean(diff(flip(Rvvb)))*ones(size(Rvvb));
dZbot = (Wall(ind(1)+1,2)-Wall(1,2))*ones(size(Zvvb));

% Output
%W = [[Rvvi;Rvvt;Rvvo;Rvvb],[-0.265906976744186;Zvvi(2:end-1);0.265906976744186;Zvvt;0.265906976744186;Zvvo(2:end-1);-0.265906976744186;Zvvb]];
W = [[Rvvi;Rvvt;Rvvo;Rvvb],[Zvvi;Zvvt';Zvvo;Zvvb']];
%dRv = [5e-3*ones(size(Rvvi)); dlh*ones(size(Rvvt));   5e-3*ones(size(Rvvo)); dlh*ones(size(Rvvb))];
dRv = [dRin; dRtop; dRout; dRbot];
%dZv = [2*0.006093023255814; dlv*ones(length(Zvvi)-2,1); 2*0.006093023255814; 12e-3*ones(size(Zvvt)); 2*0.006093023255814; dlv*ones(length(Zvvo)-2,1); 2*0.006093023255814; 12e-3*ones(size(Zvvb))];
dZv = [dZin; dZtop'; dZout; dZbot'];

%W(:,1) = circshift(W(:,1),-22);
%W(:,2) = circshift(W(:,2),-22);
%dRv    = circshift(dRv,-22);
%dZv    = circshift(dZv,-22);

end