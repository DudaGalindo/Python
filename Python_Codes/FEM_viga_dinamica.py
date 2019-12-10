''' Código para análise dinâmica de vigas'''

import numpy as np
from scipy.linalg import eig
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

    def Kg_Mg_Cg(xCC,valor_CC,Kg,Mg,Cg,ngl_tot):
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
                        if CC[j]==0:
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

        Mel = rho*A*he/420*m1 + rho*I/(30*he)*m2
        return Mel

    def Global(E,rho,A,I,cf,Fc,x,Kg,Fg,Mg,Cg,n_el,ngl_el,conec,q,alpha,beta):
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Viga.Elem(E,I,he,cf,q,x[el],x[el+1])
            Mel = VigaDin.ElemMass(rho,A,he,I)
            Cel = alpha*Kel + beta*Mel #NÃO ENTENDI PORQUÊ

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

    def deslocamento(E,rho,A,I,cf,x,xCC,valor_CC,xF,F,n_el,q,alpha,beta):
        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no)

        Fc = np.zeros(ngl_tot)
        Fc[xF[:]] = F[:]

        Kg,Fg,Mg,Cg = VigaDin.initialize(ngl_tot)
        Kg,Fg,Mg,Cg = VigaDin.Global(E,rho,A,I,cf,Fc,x,Kg,Fg,Mg,Cg,n_el,ngl_el,conec,q,alpha,beta)
        Kg,Mg,Cg = VigaDin.Kg_Mg_Cg(xCC,valor_CC,Kg,Mg,Cg,ngl_tot)

        ''' Encontrando os autovalores e autovetores '''
        lamb,psi = eig(Kg,Mg)
        w = lamb**(0.5); #w = w.real
        #psi = psi.real

        '''Tornando a matriz modal psi ortonormal à de massa'''
        MgD = (psi.T)@Mg@psi
        PSI = np.zeros(psi.shape)
        for r in range(0,ngl_tot):
            PSI[:,r] = psi[:,r]/(MgD[r,r])**(0.5)

        ''' Pré e Pós multiplicando os termos'''
        Mg = (PSI.T)@Mg@PSI
        Cg = (PSI.T)@Cg@PSI
        Kg = (PSI.T)@Kg@PSI

        '''Zerando os termos fora das diagonais (já são quase 0)'''
        for i in range(0,ngl_tot):
            for j in range(0,ngl_tot):
                if j!=i:
                    Mg[i,j] = 0
                    Kg[i,j] = 0
                    Cg[i,j] = 0

        ''' Aplicando as condições de contorno novamente ''' # se não fizer isso não da certo
        Kg,Mg,Cg = VigaDin.Kg_Mg_Cg(xCC,valor_CC,Kg,Mg,Cg,ngl_tot)

        ''' Setting for plot '''
        X = np.zeros(ngl_tot)
        x = np.linspace(0,1,ngl_tot)
        fig,ax = plt.subplots(1)
        line, = plt.plot(x, X, 'r-')
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)

        ''' Computing X '''
        for i in range(0,len(w)):
            D = -((w[i])**2)*Mg + 1j*w[i]*Cg + Kg #solve problem: D is returning as a singular matrix
            H = PSI@np.linalg.inv(D)@PSI.T
            for a in range(0,len(H[:,0])):
                for b in range(0,len(H[0,:])):
                    H[a,b] = np.linalg.norm(H[a,b])

            H = H.real
            X = H@Fg
            line.set_ydata(X.T*100)
            fig.canvas.draw_idle()
            fig.canvas.start_event_loop(.05)
            plt.show(block=False)
