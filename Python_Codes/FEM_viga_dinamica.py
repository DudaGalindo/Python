''' Código para análise dinâmica de vigas'''

import numpy as np
from .Eigen import eigV
from .FEM_viga import Viga
from .FEM import general
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class VigaDin:

    def initialize(ngl_tot):
        Mg = np.zeros([ngl_tot,ngl_tot])
        Cg = np.zeros([ngl_tot,ngl_tot])
        Kg,Fg = general.initialize(ngl_tot)
        return Kg,Fg,Mg,Cg

    def aplicacao_CC(xCC,valor_CC,Kg,Mg,Cg,ngl_tot):
        CC = np.zeros(ngl_tot)

        for i in range(0,len(xCC)):
            CC[xCC[i]] = 1
            Kg[xCC[i],xCC[i]] = 1
            Mg[xCC[i],xCC[i]] = 1
            Cg[xCC[i],xCC[i]] = 1

        for i in range(0,ngl_tot):
            if CC[i] == 1:
                for j in range(0,ngl_tot):
                    if j!=i:
                        Kg[i,j] = 0
                        Mg[i,j] = 0
                        Cg[i,j] = 0
                        if CC[j] == 0:
                            Kg[j,i] = 0
                            Mg[j,i] = 0
                            Cg[j,i] = 0
        return Kg,Mg,Cg

    def ElemMass(rho,A,he,I):
        m1 = np.array([[156,-22*he,54,13*he],\
                      [-22*he,4*he**2,-13*he,-3*he**2],\
                      [54,-13*he,156,22*he],\
                      [13*he,-3*he**2,22*he,4*he**2]])
        m2 = np.array([[36,-3*he,-36,-3*he],\
                      [-3*he,4*he**2,3*he,-he**2],\
                      [-36,3*he,36,3*he],\
                      [-3*he,-he**2,3*he,4*he**2]])

        Mel = (rho*A*he/420)*m1 + (rho*I/(30*he))*m2
        return Mel

    def Global(E,rho,A,I,cf,Fc,x,Kg,Fg,Mg,Cg,n_el,ngl_el,conec,q,alpha,beta):
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Viga.Elem(E,I,he,cf,q,x[el],x[el+1])
            Mel = VigaDin.ElemMass(rho,A,he,I)
            Cel = alpha*Kel + beta*Mel

            for i in range(0,ngl_el):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                for j in range(0,ngl_el):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]
                    Mg[ig,jg] = Mg[ig,jg] + Mel[i,j]
                    Cg[ig,jg] = Cg[ig,jg] + Cel[i,j]

        Fg = Fg + Fc
        return Kg,Fg,Mg,Cg

    def modos_vibrar(E,rho,A,I,cf,x,xCC,valor_CC,xF,F,n_el,q,alpha,beta,lastMode):
        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no)

        Fc = np.zeros(ngl_tot)
        Fc[xF[:]] = F[:]

        Kg,Fg,Mg,Cg = VigaDin.initialize(ngl_tot)
        Kg,Fg,Mg,Cg = VigaDin.Global(E,rho,A,I,cf,Fc,x,Kg,Fg,Mg,Cg,n_el,ngl_el,conec,q,alpha,beta)
        Kg,Mg,Cg = VigaDin.aplicacao_CC(xCC,valor_CC,Kg,Mg,Cg,ngl_tot)

        ''' Encontrando os autovalores e autovetores '''
        lamb,psi = eigV(Kg,Mg)
        w = lamb**(1/2)

        '''Tornando a matriz modal psi ortonormal à de massa'''
        MgD = (psi.T)@Mg@psi
        PSI = np.zeros(psi.shape)
        for r in range(0,ngl_tot):
            PSI[:,r] = psi[:,r]/(MgD[r,r])**(1/2)

        ''' Setting for plot - modos de vibrar '''

        time = np.linspace(1,1.1,1/0.01)
        modo = np.zeros(n_el+1)
        fig,ax = plt.subplots(1)
        line, = plt.plot(x, modo, 'r-')
        ax.set_xlim(0,x[n_el])
        ax.set_ylim(-1, 1)

        ''' Modos de vibração '''

        for p in range(2,lastMode):
            PSI_desloc = np.zeros(n_el+1)
            for j in range(0,n_el+1):
                PSI_desloc[j] = PSI[j*2,p]
            for t in range(0,len(time)):
                modo = PSI_desloc*np.cos(w.real[p]*time[t])
                line.set_ydata(modo*2)
                fig.canvas.draw_idle()
                fig.canvas.start_event_loop(.005)
                plt.show(block=False)

        return PSI,w,Kg,Mg,Fg,Cg


    def resposta_frequencia(PSI,w,Kg,Mg,Fg,Cg,ngl_tot,x,xCC,valor_CC,xFc):

        ''' Pré e Pós multiplicando os termos '''
        Mg = (PSI.T)@Mg@PSI
        Cg = (PSI.T)@Cg@PSI
        Kg = (PSI.T)@Kg@PSI

        ''' Zerando os valores fora das diagonais (já são quase 0) '''
        for i in range(0,ngl_tot):
            for j in range(0,ngl_tot):
                if j!=i:
                    Mg[i,j] = 0
                    Kg[i,j] = 0
                    Cg[i,j] = 0

        ''' Aplicando as condições de contorno novamente ''' # se não fizer isso não da certo
        Kg,Mg,Cg = VigaDin.aplicacao_CC(xCC,valor_CC,Kg,Mg,Cg,ngl_tot)

        n_freq = 20000
        we = np.linspace(0,n_freq,n_freq+1)
        X = np.zeros(n_freq+1)

        ''' Computing X '''
        ngl_resp = 2 #número de grau de liberdade que você deseja obter a resposta
        Hij = 0
        for w_exc in range(0,n_freq+1):

            for r in range(0,ngl_tot):
                Hij = Hij + PSI[xFc,r]*PSI[ngl_resp,r]/(-w_exc**2+1j*w_exc*Cg[r,r]+Kg[r,r])

            Hij = abs(Hij)

            for i in range(0,ngl_tot):
                X[w_exc] = X[w_exc] + Hij*Fg[i]
            Hij = 0

        plt.figure(2)
        plt.loglog(we, X, 'r-')
        plt.show()
