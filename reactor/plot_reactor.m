function plot_reactor(varargin)

%% Set up parameters
  
  % --- Define default parameters ---
  params=struct('fig',1);
  
  % --- Overwrite default parameters ---
  for j=1:length(varargin)-1 
    if ischar(varargin{j}) && any(strcmpi(fieldnames(params),varargin{j}));
      params.(lower(varargin{j}))=varargin{j+1};
    end
  end

%% Loading machine coils
W   = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_vacuum_vessel');
EF  = load('/home/fmsalvador/matlab/plasma_physics_project/reactor/reactor_coils');

%% Plotting
figure(params.fig)
grid on
hold on
xlabel('R (m)')
ylabel('Z (m)')
axis equal

% Toroidal field coils
R0 = 2.8;
B0 = 6;

% Geo params
R1 = 1.1;
R2 = 5;
d1 = 0.5;
d2 = 0.4;

[r1,z1,s1,b1] = TF_coil_shape(R1,R2,R0,B0);
[r2,z2,s2,b2] = TF_coil_shape(R1-d1,R2+d2,R0,B0);

R = [r1 r2 r1(1)];
Z = [z1 z2 z1(1)];
fill(R,Z,[1 1 1]*0.9,'edgecolor','none')
plot(r1,z1,'k','linewidth',2)
plot(r2,z2,'k','linewidth',2)

% OH Coils
for iecoils = 1:7
    plot([EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2],[EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2],'k')
end

% E, F and D Coils
for ifcoils = 8:size(EF,1)
    plot([EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2],[EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2],'r')
end

% Vacuum vessel
ind = find(isnan(W(:,1))==1);
fill([W(1:ind(1)-1,1); W(ind(1)+1:ind(2)-1,1); W(1,1)],[W(1:ind(1)-1,2); W(ind(1)+1:ind(2)-1,2); W(1,2)],[1 1 1]*0.5,'edgecolor','none')
fill([W(ind(1)+1:11,1); W(ind(2)+1:end,1); W(ind(1)+1,1)],[W(ind(1)+1:11,2); W(ind(2)+1:end,2); W(ind(1)+1,2)],[1 1 1]*0.8,'edgecolor','none')
plot(W(:,1),W(:,2),'k','linewidth',1)

end