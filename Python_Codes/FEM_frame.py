import numpy as np
import math
from .FEM import general

class Frame:
    def Elem(A,E,beta,he,I):
        c = math.cos(beta); s = math.sin(beta)
        a1 = 4*E*I/he; a2 = 2*E*I/he; a3 = A*E*c**2/he + 12*E*I*s**2/he**3;
        a4 = A*E*c*s/he - 12*E*I*c*s/he**3; a5 = -6*E*I*s/he**2;
        a6 = A*E*s**2/he + 12*E*I*c**2/he**3; a7 = 6*E*I*c/he**2
        Kel = np.array([[a3,a4,a5,-a3,-a4,a5],[a4,a6,a7,-a4,-a6,a7],\
                        [a5,a7,a1,-a5,-a7,a2],[-a3,-a4,-a5,a3,a4,-a5],\
                        [-a4,-a6,-a7,a4,a6,-a7],[a5,a7,a2,-a5,-a7,a1]])

        return Kel

    def Global(conec,conec_el,n_el,coord_no,ngl_el,A,E,I,Kg):
        for el in range(0,n_el):
            No1 = int(conec_el[el,0])
            No2 = int(conec_el[el,1])
            L = math.sqrt((coord_no[No2-1,1] - coord_no[No1-1,1])**2+(coord_no[No2-1,0] - coord_no[No1-1,0])**2)
            he = L
            if (coord_no[No2-1,0] - coord_no[No1-1,0])==0:
                beta = 2*math.atan(1)
            else: beta = math.atan((coord_no[No2-1,1] - coord_no[No1-1,1])/(coord_no[No2-1,0] - coord_no[No1-1,0]))
            Kel = Frame.Elem(A,E,beta,he,I)
            for i in range(0,ngl_el):
                ig = int(conec[el,i] - 1)
                for j in range(0,ngl_el):
                    jg = int(conec[el,j] - 1)
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]
        return Kg

    def deslocamento(coord_no,n_nos_el,F,xF,xCC,valor_CC,n_el,A,E,I):
        n_nos_tot = general.n_nosTOTAL(n_nos_el,n_el)
        ngl_no = 3*np.ones(n_nos_tot)
        ngl_no = ngl_no.astype(int)
        ngl_el, ngl_tot = general.ngl(ngl_no[0],n_nos_el[0],n_nos_tot)

        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no[0])
        conec_el = general.conect(1*n_nos_tot,n_el,2,1)

        Kg,Fg = general.initialize(ngl_tot)

        for i in range(len(xF)):
            Fg[xF[i]] = F[i]

        Kg = Frame.Global(conec,conec_el,n_el,coord_no,ngl_el,A,E,I,Kg)
        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg = np.linalg.inv(Kg)
        u = np.matmul(Kg,Fg.T)
        return u
