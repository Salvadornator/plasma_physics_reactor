function gradinho(r0,z0,psi0,J0,I_Fcoils)

% vessel walls
paredes_r = [0.4000 0.4000 0.8450 0.8450 0.4000];
paredes_z = [-0.2600 0.2600 0.2600 -0.2600 -0.2600];

% mask for arrays inside vessel
r0_vessel_mask=(min(paredes_r)<=r0) & (max(paredes_r)>=r0);
z0_vessel_mask=(min(paredes_z)<=z0) & (max(paredes_z)>=z0);
r0_vessel = r0(r0_vessel_mask);
z0_vessel = z0(z0_vessel_mask);

% create new r and z inside vessel
tamanho = 65;
r_vessel = linspace(min(paredes_r), max(paredes_r), tamanho);
z_vessel = linspace(min(paredes_z), max(paredes_z), tamanho);

% interpolate psi and J values inside vessel
psi_vessel = interp2(r0, z0, psi0, r_vessel, z_vessel,'cubic');
J_vessel = interp2(r0, z0, J0, r_vessel, z_vessel, 'cubic');
size(psi_vessel)

% create D matrix. f is the poisson equation source, and also returns dr and dz
[D, f, dr, dz] = create_D(r_vessel, z_vessel, psi_vessel, J_vessel);


% xs, Dsinv, ps, qs = bicicleta(D, f)
% psi = back_sub(xs, Dsinv, ps, qs, psi_vessel)
% sys.exit()
% # np.save("grid_inner_vessel_65",[r_vessel,z_vessel,[r0_vessel_min,r0_vessel_max],[z0_vessel_min,z0_vessel_max],dr,dz])

% J_vessel_rec = remake_J(r_vessel, z_vessel, psi)

% r0_vessel_min = np.where(r0_vessel_mask)[0].min()
% r0_vessel_max = np.where(r0_vessel_mask)[0].max()
% z0_vessel_min = np.where(z0_vessel_mask)[0].min()
% z0_vessel_max = np.where(z0_vessel_mask)[0].max()


% psi_reconstructed = redraw_field(r0_vessel, z0_vessel, r_vessel, z_vessel, psi0, psi, [
%                                  r0_vessel_min, r0_vessel_max], [z0_vessel_min, z0_vessel_max])

% J_reconstructed = redraw_field(r0_vessel, z0_vessel, r_vessel[1:-1], z_vessel[1:-1], J0, J_vessel_rec, [
%     r0_vessel_min, r0_vessel_max], [z0_vessel_min, z0_vessel_max])

% print("Tempo para calcular: {:.3f} ms".format(
%     1000 * (time.time() - time0)))

% ploti = 1
% if ploti == 1:
%     fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6, 6))
%     # fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(14, 6))
%     im0 = ax[0, 0].contour(r0, z0, psi0, cmap="plasma",
%                            levels=np.linspace(0, 0.23, 30))
%     ax[0, 0].set_title("Original")
%     fig.colorbar(im0, ax=ax[0, 0])
%     im1 = ax[0, 1].contour(r0, z0, psi_reconstructed, cmap="plasma",
%                            levels=np.linspace(psi.min(), psi.max(), 30))
%     # im1 = ax[0, 1].imshow(psi_reconstructed, cmap="plasma",
%     #                        )        
%     # ax[0, 1].set_title("Reconstruído")
%     fig.colorbar(im1, ax=ax[0, 1])
%     # im2 = ax[2].imshow(J_source)
%     im2 = ax[1, 0].contour(r_vessel, z_vessel, J_source,
%                            cmap="plasma", levels=np.linspace(0, 3e6, 20))
%     # im2 = ax[1, 0].imshow(J_source)
%     fig.colorbar(im2, ax=ax[1, 0])
%     ax[1, 0].set_title("Densidade \nde Corrente")
%     ax[1, 0].set_xlim(r_vessel[0], r_vessel[-1])
%     ax[1, 0].set_ylim(z_vessel[0], z_vessel[-1])
%     # ,levels=np.linspace(0,3e6,20))
%     im3 = ax[1, 1].contour(
%         r_vessel[1:-1], z_vessel[1:-1], J_vessel_rec, cmap="plasma", levels=np.linspace(0, 3e6, 20))
%     # im3 = ax[1, 1].imshow(J_reconstructed)
%     ax[1, 1].set_title("Densidade de Corrente \n Reconstruída")
%     ax[1, 1].set_xlim(r_vessel[0], r_vessel[-1])
%     ax[1, 1].set_ylim(z_vessel[0], z_vessel[-1])
%     fig.colorbar(im3, ax=ax[1, 1])
%     ax[0, 0].plot(paredes_r, paredes_z, color="k")
%     ax[0, 1].plot(paredes_r, paredes_z, color="k")
%     ax[1, 0].plot(paredes_r, paredes_z, color="k")
%     ax[1, 1].plot(paredes_r, paredes_z, color="k")
%     # ax[0, 0].set_aspect('equal', adjustable='box')
%     # ax[0, 1].set_aspect('equal', adjustable='box')
%     # ax[1, 0].set_aspect('equal', adjustable='box')
%     # ax[1, 1].set_aspect('equal', adjustable='box')
%     fig.tight_layout()

%     # plt.savefig(config)
% else:
%     om.io.savemat("rec_fields",{"J_rec":J_reconstructed,"psi_rec":psi_reconstructed})
%     rec_fields=om.io.loadmat("rec_fields.mat")
%     J_vessel_rec=rec_fields["J_rec"]
%     plt.contour(
%         r0, z0, J_vessel_rec, cmap="plasma", levels=np.linspace(0, 3e6, 20))


% # print(D)
% # print(f)
% return r_vessel, z_vessel, psi_vessel
