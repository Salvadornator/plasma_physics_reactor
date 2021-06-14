"""
Author: C. Amador 2019/06
Finds X and O points given a plasma psi profile
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import measure


def XOXO(r, z, psi):
    """procura-se pontos em que primeira derivada de R e Z são 0 ao mesmo tempo
    esses pontos são analisado na segunda derivada, ver se sinal é igual (multiplica e vê resultado)
    ponto O sinal é igual
    ponto X é diferente
    """
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    tamanho = r.size
    psi1_r = np.gradient(psi, dr, axis=1)
    psi1_z = np.gradient(psi, dz, axis=0)
    psi2_r = np.gradient(psi1_r, dr, axis=1)
    psi2_z = np.gradient(psi1_z, dz, axis=0)
    psi2_zr = np.gradient(psi1_z, dr, axis=1)
    B_ini = np.sqrt(psi1_r**2 + psi1_z**2)
    # ponto O
    ploti = 0
    B = B_ini[int(tamanho * 0.25):int(tamanho * 0.75),
              int(tamanho * 0.25):int(tamanho * 0.75)]
    y, x = np.where(B == B.min())
    rE = r[int(tamanho * 0.25):int(tamanho * 0.75)][x]
    zE = z[int(tamanho * 0.25):int(tamanho * 0.75)][y]
    if ploti == 1:
        plt.figure()
        plt.imshow(B_ini)
        plt.plot(x + tamanho * 0.25, y + tamanho * 0.25, "ko", lw=5)
    fE = psi[y, x]

    # ponto X
    rS = []
    zS = []
    fS = []
    D = (psi2_r * psi2_z) - (psi2_zr**2)
    mask=np.where(D[1:-1,1:-1] < 0)
    for m in range(mask[0].size):
        j=mask[0][m]+1
        i=mask[1][m]+1
        B_min=np.min(B_ini[j-1:j+2,i-1:i+2])
        if B_min==B_ini[j,i]:
            rX = r[i]
            zX = z[j]
            rS.append(rX)
            zS.append(zX)
            fS.append(psi[j,i])
            if ploti == 1:
                plt.plot(i,j, "rx", lw=5)
    
    psi_borda=psi.copy()
    psi_borda[1:-1,1:-1]=0
    contours = measure.find_contours(psi, psi_borda.max())
    psi_lcfs=psi_borda.max()

    if len(fS)!=0:
        psiX=0
        for ffS in fS:
            contours = measure.find_contours(psi, ffS*0.995)
            for n, contour in enumerate(contours):
                if np.all(contour[0]==contour[-1]):
                    psiX=ffS*0.995
        if psiX>psi_lcfs:
            psi_lcfs=psiX[0]
  
    psin = (psi - fE) / (psi_lcfs - fE)

    if ploti==1:
        # # Display the image and plot all contours found
        # fig, ax = plt.subplots()
        # ax.imshow(psi, interpolation='nearest', cmap=plt.cm.gray,origin='upper')

        contours = measure.find_contours(psi, psi_lcfs)

        for n, contour in enumerate(contours):
            plt.plot(contour[:, 1], contour[:, 0], linewidth=2)

        # ax.axis('image')
        plt.xticks([])
        plt.yticks([])
        plt.show()


    return 0, 0, 0, 0, 0, 0, 0

if __name__ == '__main__':
    import grad_shafranov_slim as gs
    # "DN", "limited", "LSN1", "LSN2", "SF", "USN"
    r, z, x = gs.toka_aki("SF")
    plt.show()

"""
psi normalizado tem que analisar a lcfs. 
interpola e procura o maior valor da na borda do vaso, e veja se tem linha fechada nele. 
calcula ponto X e vê se o valor de psi é maior que ele. O ponto X precisa estar numa linha fechada para contar.
"""