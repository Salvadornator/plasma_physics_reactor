%salvar config nova como equil_profiles/run_equil.mat
% ela deve ter: r,z,psi,J,I_Fcoils
tic;
system('/opt/anaconda5/bin/python3.6 grad_shafranov.py');
fields=load('rec_fields');
J_rec=fields.J_rec;
psi_rec=fields.psi_rec;
rE=fields.rE;
zE=fields.zE;
fE=fields.fE;
rS=fields.rS;
zS=fields.zS;
fS=fields.fS;
psin=fields.psin;


toc
% equilibrium=load('equil_profiles/run_equil.mat');
% r=equilibrium.r;
% z=equilibrium.z;
% psi=equilibrium.psi;
% J=equilibrium.J;
% I_Fcoils=equilibrium.I_Fcoils;
% gradinho(r,z,psi,J,I_Fcoils);
