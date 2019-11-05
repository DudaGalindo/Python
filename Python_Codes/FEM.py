import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from scipy.integrate import quad

class general:
    def n_nosTOTAL(n_nos_el,n_el):
        n_nos_total = int(sum(n_nos_el) - (n_el-1))
        return n_nos_total

    def ngl(ngl_no,n_nos_el,n_nos_tot):
        ngl_el = ngl_no*n_nos_el
        ngl_el = ngl_el.astype(int)
        ngl_tot = ngl_no*n_nos_tot #número de GL da estrutura
        ngl_tot = ngl_tot.astype(int)
        return ngl_el,ngl_tot

    def initialize(ngl_tot):
        Kg = np.zeros([ngl_tot,ngl_tot])
        Fg = np.zeros(ngl_tot)
        return Kg, Fg

    def conect(ngl_tot,n_el,ngl_el,ngl_no):
        conec = np.zeros((n_el,ngl_el));
        for el in range(0,n_el):
            for i in range(0,ngl_el):
                conec[el,i] = i+ngl_no*el+1
        conec = conec.astype(int)
        return conec

    def Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg): ##acho que tem como otimizar
        CC = np.zeros(ngl_tot)

        for i in range(0,len(xCC)):
            CC[xCC[i]] = 1
            Fg[xCC[i]] = valor_CC[i]
            Kg[xCC[i],xCC[i]] = 1

        for i in range(0,ngl_tot):
            if CC[i] == 1:
                for j in range(0,ngl_tot):
                    if j!=i:
                        Kg[i,j] = 0
                        if CC[j]==0:
                            Fg[j] = Fg[j] - Kg[j,i]*Fg[i]
                            Kg[j,i] = 0
        return Kg,Fg
