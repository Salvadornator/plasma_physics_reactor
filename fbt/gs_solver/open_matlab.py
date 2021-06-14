import numpy as np
from scipy import io
import matplotlib.pyplot as plt


def current(r, z, I_Fcoils):
    folder = "equil_reconst/gs/"
    r0, z0, dR, dZ, Nturn = np.loadtxt(
        folder + "tcabr_coils", skiprows=8, usecols=range(5), unpack=True)
    J = I_Fcoils / (dR * dZ)
    J_Fcoils = np.zeros((r.size, z.size))
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    r_delta = np.array(dR / (2 * dr), dtype=int)
    z_delta = np.array(dZ / (2 * dz), dtype=int)

    for coil in range(J.size):
        r_pos = abs(r0[coil] - r).argmin()
        z_pos = abs(z0[coil] - z).argmin()
        J_Fcoils[z_pos - z_delta[coil]:z_pos + z_delta[coil],
                 r_pos - r_delta[coil]:r_pos + r_delta[coil]] = J[coil]
    # plt.title("Control Coils")
    # plt.contourf(r,z,J_Fcoils)
    # plt.show()
    return J_Fcoils


def open_equilibrium(name):
    folder = ""
    equilibrium = io.loadmat(folder + name)
    r = equilibrium['r'].flatten()
    z = equilibrium['z'].flatten()
    psi = equilibrium['psi']
    J = equilibrium['J']
    # plt.contour(r,z,psi)
    # plt.show()
    return r, z, psi, J

if __name__ == '__main__':
    name = "limited.mat"
    r, z, psi, J, I_Fcoils = open_equilibrium(name)
    J_Fcoils = current(r, z, I_Fcoils)
    # print(psi)
    fig=plt.figure()
    plot1=fig.add_subplot(1,2,1)
    plot1.contourf(r,z,J+J_Fcoils)
    # plot1.colorbar()
    plot2=fig.add_subplot(1,2,2)
    plot2.contourf(r,z,psi)
    plt.show()
