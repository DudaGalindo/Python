import numpy as np
import math
from .FEM import general

class Barra:

    def Elem(n_nos_el,integralAdx,E,he,f):

        if n_nos_el == 3:
            k = 1./3*np.array([[7., -8., 1.],[-8., 16., -8.],[1., -8., 7.]])
            fel = np.array([1./6, 2./3, 1./6])

        if n_nos_el == 2:
            k = np.array([[1., -1.],[-1., 1.]])
            fel = np.array([1./2, 1./2])

        Kel = (E/he**2)*integralAdx*k    #quad(A,x1,x2)
        Fel = f*he*fel
        return Kel,Fel

    def Global(n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,Kg,Fg):

        ##conec = np.array([[1,2],[2,3]])
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Barra.Elem(n_nos_el[el],integralAdx(he,el),E[el],he,f[el])
            for i in range(0,n_nos_el[el]):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                for j in range(0,n_nos_el[el]):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]

        Fg = Fg + Fc
        return Kg, Fg

    def deslocamento(n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC):
        n_nos_tot = general.n_nosTOTAL(n_nos_el,n_el) #número total de nós da estrutura
        ngl_no = 1*np.ones(n_el)
        ngl_no = ngl_no.astype(int)
        #ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el,n_nos_tot)
        Kg,Fg = general.initialize(n_nos_tot)
        Kg,Fg = Barra.Global(n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,Kg,Fg)
        Kg,Fg = general.Kg_Fg(n_nos_tot,xCC,valor_CC,Kg,Fg)
        Kg_inv = np.linalg.inv(Kg)
        Fg = (Fg.T)
        u = np.matmul(Kg_inv,Fg)
        return u
