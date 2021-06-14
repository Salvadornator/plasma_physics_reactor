"""
Parede externa:
0.395<R<0.850
-0.272<Z<0.272
[0.3950,0.3950,0.8500,0.8500,0.3950]
[-0.2720, 0.2720, 0.2720, -0.2720, -0.2720]
Parede interna:
0.400<R<0.845
-0.260<Z<0.260
[0.4000,0.4000,0.8450,0.8450,0.4000]
[-0.2600,0.2600,0.2600,-0.2600,-0.2600]
"""
import sys
import os
import copy
import time
import numpy as np
from scipy import constants, interpolate
import matplotlib.pyplot as plt
import matplotlib
import open_matlab as om
import xoxo

np.set_printoptions(linewidth=160)  # ,threshold=np.inf)

mu0 = constants.value('magn. constant')


def create_f(r, psi_vessel, J_source, dz, a, b):
    f = 2 * np.pi * mu0 * r * J_source * (dz**2)

    # aplicando condições de contorno
    f[1, :] += psi_vessel[0, :]
    f[-2, :] += psi_vessel[-1, :]
    f[:, 1] += b[1] * psi_vessel[:, 0]
    f[:, -2] += a[-2] * psi_vessel[:, -1]

    f = f[1:-1, 1:-1]

    return f


def create_D(r_vessel, z_vessel, psi_vessel, J_source, tamanho):
    """
    matriz
        | D  -I  0 |
    A = |-I   D -I |
        | 0  -I  D |
    """

    dr = r_vessel[1] - r_vessel[0]
    dz = z_vessel[1] - z_vessel[0]
    dz_dr2 = (dz / dr)**2

    # For central difference:
    a = dz_dr2 * r_vessel / (r_vessel + (dr / 2))
    b = dz_dr2 * r_vessel / (r_vessel - (dr / 2))
    # For forward difference:
    # a = dz_dr2 * (1 - (dr / r[1:-1]))
    # b = dz_dr2 * np.ones(r[1:-1].size)
    e = - 2 - a - b

    f = create_f(r_vessel, psi_vessel, J_source, dz, a, b)

    a = a[1:-1]
    b = b[1:-1]
    e = e[1:-1]

    D = np.zeros_like(f)

    for i in range(D.shape[0]):
        D[i, i] = - e[i]
        if 0 < i < D.shape[0] - 1:
            D[i, i - 1] = - b[i]
            D[i, i + 1] = - a[i]
        elif i == 0:
            D[i, i + 1] = - a[i]
        elif i == D.shape[0] - 1:
            D[i, i - 1] = - b[i]

    np.save("BCR_{}_b".format(tamanho), [D, dr, dz, a, b])

    return D, f, dr, dz


def eval_Ds(D, step, I):
    pow2 = (2**step)
    i = np.arange(1, (2**step) + 1)
    theta = (i - 0.5) * np.pi / pow2
    lambida_i = 2 * np.cos(theta)
    a = np.zeros_like(D) + 1
    b = np.zeros_like(D)
    for s in range(lambida_i.size):
        termo = (D - lambida_i[s] * I)
        a = np.dot(a, termo)
        b += (((-1)**s) * np.linalg.inv(termo) / pow2) * np.sin(theta[s])
    return a, b


def bicicleta(D, f, tamanho):
    """
    equação: prod (1->ks) [A-lambda_i^s*I] * y = prod (1->ls) [A-mu_i^s*I] * q_js^s 
    Case II from paper R. A. Sweet https://www.jstor.org/stable/2156489
    """
    ta_ocupado = True
    if os.path.isfile("BCR_{}_c.npy".format(tamanho)):
        try:
            Ds, Dsinv, steps = np.load("BCR_{}_c.npy".format(tamanho),allow_pickle=True)
            ta_ocupado = False
        except ValueError:
            ta_ocupado=True

    I = np.identity(D.shape[0])
    if ta_ocupado:
        steps = 1
        s = D.shape[0]
        while s != 1:
            s = int((s - 1) / 2)
            steps += 1
        Ds = np.zeros((steps, D.shape[0], D.shape[1]))
        Dsinv = np.zeros_like(Ds)
        Ds[0] = D
        Dsinv[0] = np.linalg.inv(Ds[0])
    qs = [f.copy()]
    ps = [np.zeros(f.shape)]

    for s in range(1, steps):
        pp = np.zeros((int(qs[-1].shape[0] / 2), qs[-1].shape[1]))
        qq = np.zeros((int(qs[-1].shape[0] / 2), qs[-1].shape[1]))

        # j é o índice do vetor antigo. pos é a posição no novo vetor
        for j in range(1, qs[-1].shape[0], 2):
            pos = int((j - 1) / 2)
            pp[pos] = ps[-1][j] + \
                np.dot(Dsinv[s - 1], (qs[-1][j] +
                                      ps[-1][j - 1] + ps[-1][j + 1]))
            qq[pos] = qs[-1][j - 1] + qs[-1][j + 1] + 2 * pp[pos]
            # print(pos)
        ps.append(pp)
        qs.append(qq)
        if ta_ocupado:
            a, b = eval_Ds(D, s, I)
            Ds[s], Dsinv[s] = a, b

    zs = qs[-1].flatten()
    vs = ps[-1].flatten() + np.dot(Dsinv[-1], zs)
    if ta_ocupado:
        np.save("BCR_{}_c".format(tamanho), [
                Ds, Dsinv, steps])
    return vs, Dsinv, ps, qs


def back_sub(xs, Dsinv, ps, qs, boundary):
    Xs = np.zeros_like(boundary)
    Xs[-1, :] = boundary[-1, :]
    Xs[0, :] = boundary[0, :]
    Xs[1:-1, 0] = boundary[1:-1, 0]
    Xs[1:-1, -1] = boundary[1:-1, -1]
    middle = int(ps[0].shape[0] / 2)
    centros = [middle + 1]
    Xs[centros[0], 1:-1] = xs
    teste = np.zeros(boundary.shape[0], dtype=int)

    teste[centros[0]] = centros[0]
    teste[-1] = len(teste) - 1
    for step in range(len(Dsinv) - 2, -1, -1):
        ps_atual = ps[step]
        qs_atual = qs[step]
        # print("\n step by step")
        # print(step, ps_atual.shape)
        tamanho_atual = ps_atual.shape[0]
        passo = 2**step
        centros2 = copy.copy(centros)
        pontos = []
        for centro in centros:
            for pos in range(centro - passo, centro + passo + 1, passo):
                if pos not in centros2:
                    centros2.append(pos)
                    if pos not in pontos:
                        pontos.append(pos)
        pontos.sort()
        centros2.sort()
        count = 0
        for pos in pontos:
            # print(pos, count * 2, step)
            if count == 0:
                Xs[pos, 1:-1] = ps_atual[2 * count, :] + np.dot(
                    Dsinv[step], (qs_atual[2 * count, :] + Xs[pos + passo, 1:-1]))

            elif count == int(tamanho_atual / 2):
                Xs[pos, 1:-1] = ps_atual[2 * count, :] + np.dot(
                    Dsinv[step], (qs_atual[2 * count, :] + Xs[pos - passo, 1:-1]))

            else:
                Xs[pos, 1:-1] = ps_atual[2 * count, :] + np.dot(
                    Dsinv[step], (qs_atual[2 * count, :] + Xs[pos + passo, 1:-1] + Xs[pos - passo, 1:-1]))
            count += 1
            teste[pos] = pos
        # print(pontos)
        # print(centros2)
        # print("passo",passo)
        # print(teste)
        # print(Xs[middle+1],Xs[middle+2])
        # print(Xs)
        # input()
        centros = copy.copy(centros2)
    # return negative because the whole poisson equation was multiplied by -1
    return Xs


def remake_J(r, z, psi):
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    psi1_r = np.gradient(psi, dr, axis=1)
    psi2_r = np.gradient(psi1_r, dr, axis=1)
    psi1_z = np.gradient(psi, dz, axis=0)
    psi2_z = np.gradient(psi1_z, dz, axis=0)
    dpsi2_z = psi2_z
    # dpsi2_r = (psi2_r / r[1:-1])
    dpsi2_r = psi2_r

    dpsi1_r = psi1_r

    parte_r = dpsi2_r - (dpsi1_r / r)
    lhs = dpsi2_z + parte_r
    J = -lhs / (mu0 * r * 2 * np.pi)
    return J


def redraw_field(r0_vessel, z0_vessel, r_vessel, z_vessel, field0, field_vessel, limits_r, limits_z):
    """
    Paste the new evaluated field into the old one grid
    """
    field = field0.copy()
    func_field = interpolate.interp2d(
        r_vessel, z_vessel, field_vessel, kind="cubic")
    field[limits_z[0]:limits_z[1] + 1, limits_r[0]
        :limits_r[1] + 1] = func_field(r0_vessel, z0_vessel)
    return field



def multi_conf():
    configs = "DN", "limited", "LSN1", "LSN2", "SF", "USN"
    for config in configs:
        print(config)
        toka_aki(config)


def toka_aki(config="run_equil"):
    time0 = time.time()
    name = config + ".mat"
    r0, z0, psi0, J0 = om.open_equilibrium(name)
    tamanho = r0.size
    ta_ocupado=True
    if os.path.isfile("BCR_{}_b.npy".format(tamanho)):
        try:
            D, dr, dz, a, b = np.load("BCR_{}_b.npy".format(tamanho),allow_pickle=True)
            f = create_f(r0, psi0, J0, dz, a, b)
            ta_ocupado=False
        except ValueError:
            print("\n Avisa o Cassio que o npy precisa ser arrumado!\n Isso não afeta o cálculo atual...\n")
            ta_ocupado=True
    if ta_ocupado:
        D, f, dr, dz = create_D(
            r0, z0, psi0, J0, tamanho)

    xs, Dsinv, ps, qs = bicicleta(D, f, tamanho)

    psi_reconstructed = back_sub(xs, Dsinv, ps, qs, psi0)

    J_vessel_rec = remake_J(r0, z0, psi_reconstructed)

    J_reconstructed = J_vessel_rec

    # J_reconstructed = redraw_field(r0_vessel, z0_vessel, r_vessel[1:-1], z_vessel[1:-1], J0, J_vessel_rec, [
    #     r0_vessel_min, r0_vessel_max], [z0_vessel_min, z0_vessel_max])

    # sys.exit()

    #rE, zE, fE, rS, zS, fS, psin = xoxo.XOXO(r0, z0, psi_reconstructed)

    ploti = 0
    if ploti == 1:
        print("Tempo para calcular: {:.3f} ms".format(
            1000 * (time.time() - time0)))
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6, 6))
        # fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(14, 6))
        im0 = ax[0, 0].contour(r0, z0, psi_0, cmap="plasma",
                               levels=np.linspace(psi_0.min(), psi_0.max(), 30))
        ax[0, 0].set_title("Psi Original")
        fig.colorbar(im0, ax=ax[0, 0])
        im1 = ax[0, 1].contour(r0, z0, psi_reconstructed, cmap="plasma",
                               levels=np.linspace(psi_0.min(), psi_0.max(), 30))
        ax[0, 1].plot(pontos[0], pontos[1], "x", lw=5)
        # im1 = ax[0, 1].imshow(psi_reconstructed, cmap="plasma",
        #                        )
        ax[0, 1].set_title("Psi Reconstruído")
        fig.colorbar(im1, ax=ax[0, 1])
        # im2 = ax[2].imshow(J_source)
        im2 = ax[1, 0].contour(r0, z0, J0,
                               cmap="plasma", levels=np.linspace(0, 3e6, 20))
        # im2 = ax[1, 0].imshow(J_source)
        fig.colorbar(im2, ax=ax[1, 0])
        ax[1, 0].set_title("Densidade \nde Corrente")
        ax[1, 0].set_xlim(r_vessel[0], r_vessel[-1])
        ax[1, 0].set_ylim(z_vessel[0], z_vessel[-1])
        # ,levels=np.linspace(0,3e6,20))
        im3 = ax[1, 1].contour(
            r0, z0, J_reconstructed, cmap="plasma", levels=np.linspace(0, 3e6, 20))
        # im3 = ax[1, 1].imshow(J_reconstructed)
        ax[1, 1].set_title("Densidade de Corrente \n Reconstruída")
        ax[1, 1].set_xlim(r_vessel[0], r_vessel[-1])
        ax[1, 1].set_ylim(z_vessel[0], z_vessel[-1])
        fig.colorbar(im3, ax=ax[1, 1])
        # ax[0, 0].plot(paredes_r, paredes_z, color="k")
        # ax[0, 1].plot(paredes_r, paredes_z, color="k")
        ax[1, 0].plot(paredes_r, paredes_z, color="k")
        ax[1, 1].plot(paredes_r, paredes_z, color="k")
        # ax[0, 0].set_aspect('equal', adjustable='box')
        # ax[0, 1].set_aspect('equal', adjustable='box')
        # ax[1, 0].set_aspect('equal', adjustable='box')
        # ax[1, 1].set_aspect('equal', adjustable='box')
        fig.tight_layout()

        # plt.savefig(config)
        plt.show()
    else:
        om.io.savemat(
            "rec_fields.mat", {"J_rec": J_reconstructed, "psi_rec": psi_reconstructed})
        # rec_fields = om.io.loadmat("rec_fields.mat")

    return r0, z0, psi_reconstructed


if __name__ == '__main__':
    r, z, x = toka_aki()
