''' Código para análise dinâmica de vigas'''

import numpy as np
import scipy as scp
from FEM import general
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class VigaDin:
    def initialize():
        Mg = np.zeros([ngl_tot,ngl_tot])
        Cg = np.zeros([ngl_tot,ngl_tot])
        Kg,Fg = general.initialize(ngl_tot)
        return Kg,Fg,Mg,Cg

    def Kg_Fg_Mg(xCC,valor_CC,Kg,Fg,Mg):
        CC = np.zeros(ngl_tot)

        for i in range(0,len(xCC)):
            CC[xCC[i]] = 1
            Fg[xCC[i]] = valor_CC[i]
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
                            Fg[j] = Fg[j] - Kg[j,i]*Fg[i]
                            Kg[j,i] = 0
                            Mg[j,i] = 0
                            Cg[j,i] = 0
        return Kg,Fg,Mg

    def ElemMass(ro,A,he,I):
        m1 = np.array([[156,-22*he,54,13*he],\
                      [-22*he,4*he**2,-13*he,-3*he**2],\
                      [54,-13*he,156,22*he],\
                      [13*he,-3*he**2,22*he,4*he**2]])
        m2 = np.array([[36,-3*he,-36,-3*he],\
                      [-3*he,4*he**2,3*he,-he**2],\
                      [-36,3*he,36,3*he],\
                      [3*he,-he**2,3*he,4*he**2]])
        Mel = ro*A*he/420*m1 + ro*I/(30*he)*m2
        return Mel

    def Global(E,ro,I,cf,Fc,x,Kg,Fg,Mg,Cg,n_el,ngl_el,conec,q,alpha,beta):
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Viga.Elem(E,I,he,cf,q,x[el],x[el+1])
            Mel = VigaDin.ElemMass(ro,A,he,I)
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

    def deslocamento(E,I,cf,x,xCC,valor_CC,xF,F,n_el,q,alpha,beta):
        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no)

        Fc = np.zeros(ngl_tot)
        Fc[xF[:]] = F[:]

        Kg,Fg,Mg,Cg = VigaDin.initialize(ngl_tot)
        Kg,Fg,Mg,Cg = Viga.Global(E,I,cf,Fc,x,Kg,Fg,Mg,n_el,ngl_el,conec,q,alpha,beta)
        Kg,Fg,Mg,Cg = VigaDin.Kg_Fg_Mg(xCC,valor_CC,Kg,Fg,Mg,Cg)

        ''' Encontrando os autovalores e autovetores '''
        lamb,phi = scp.linalg.eig(Kg*inv(Mg))
        w = lamb**(1/2)

        '''Tornando a matriz modal psi ortonormal à de massa'''
        PSI = np.zeros(psi.shape)
        for r in range(psi[:,0]):
            PSI[:,r] = psi[:,r]/sqrt(Mg[r,r])

        ''' Setting for plot '''
        X = np.zeros(len(Fg))
        fig,ax = plt.subplots(1)
        line, = plt.plot(x, X, 'r-')
        ax.set_xlim(0, x[n_el-1])
        ax.set_ylim(-1, 1)

        for i in range(len(w)):
            D = -w[i]**2*PSI.T*Mg*PSI + 1j*w[i]*PSI.T*C*PSI + PSI.T*Kg*PSI
            H = PSI*inv(D)*PSI.T
            X = H*Fg

            line.set_ydata(X)
            fig.canvas.draw_idle()
            fig.canvas.start_event_loop(.05)
            plt.show(block=False)

        #lamb = sp.symbols('lambda')
        #cp = sp.det(Kg-lamb*Mg)
        #eigs = sp.solveset(cp,lamb)
        #eigs_v2 = sp.solveset(cp.rewrite(sp.exp).simplify(), lam)
