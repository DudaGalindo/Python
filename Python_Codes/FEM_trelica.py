import numpy as np
import math
from .FEM import general


class Trelica:
    def conect(n_el,ngl_tot,ngl_el,coord_no,conec_el): ##TENTAR OTIMIZAR
        conec = np.zeros([n_el,ngl_el])
        vetor = np.linspace(1,ngl_tot,ngl_tot)
        for el in range(0,n_el):
            i = 0
            for No in range(0,ngl_el):
                i = np.sign(No-i); c = np.sign(abs(2-No)*No)
                conec[el,No] = vetor[(conec_el[el,i]-1)*2 + c]
                i = i+1
        conec = conec.astype(int)
        return conec

    def Elem(A,E,beta,he):
        c = math.cos(beta); s = math.sin(beta)
        k = np.array([[c*c, c*s, -c*c, -c*s],[c*s, s*s, -c*s, -s*s],[-c*c, -c*s, c*c, c*s],[-c*s, -s*s, c*s, s*s]])
        Kel = A*E/he*k
        return Kel

    def Melem(rho,A,L,beta):
        c = math.cos(beta); s = math.sin(beta)
        #matriz de massa consistente:
        m = np.array([[2*c*c,2*c*s, c*c, c*s],[2*c*s, 2*s*s, c*s, s*s],[c*c, c*s, 2*c*c, 2*c*s], [c*s, s*s, 2*c*s, 2*s*s]])
        Mel = rho*A*L/6*m
        return Mel

    def Global(conec,conec_el,n_el,coord_no,ngl_el,A,E,Kg):
        for el in range(0,n_el):
            No1 = int(conec_el[el,0])
            No2 = int(conec_el[el,1])
            L = math.sqrt((coord_no[No2-1,1] - coord_no[No1-1,1])**2+(coord_no[No2-1,0] - coord_no[No1-1,0])**2)
            he = L
            if (coord_no[No2-1,0] - coord_no[No1-1,0])==0:
                beta = 2*math.atan(1)
            else: beta = math.atan((coord_no[No2-1,1] - coord_no[No1-1,1])/(coord_no[No2-1,0] - coord_no[No1-1,0]))
            Kel = Trelica.Elem(A[el],E[el],beta,he)
            for i in range(0,ngl_el):
                ig = conec[el,i] - 1
                for j in range(0,ngl_el):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]
        return Kg

    def deslocamento_tensao(coord_no,conec_el,n_nos_tot,n_nos_el,F,xF,xCC,valor_CC,n_el,A,E):
        ngl_no = 2*np.ones(n_nos_tot)
        ngl_no = ngl_no.astype(int)
        ngl_el, ngl_tot = general.ngl(ngl_no[0],n_nos_el[0],n_nos_tot)

        conec = Trelica.conect(n_el,ngl_tot,ngl_el,coord_no,conec_el)

        # Deslocamento:
        Kg,Fg = general.initialize(ngl_tot)

        for i in range(len(xF)):
            Fg[xF[i]] = F[i]

        Kg = Trelica.Global(conec,conec_el,n_el,coord_no,ngl_el,A,E,Kg)
        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg = np.linalg.inv(Kg)
        u = np.matmul(Kg,Fg.T)

        # Tensões:
        u_el = np.zeros((ngl_el,1))
        Fel = np.zeros(ngl_el)
        Stress = np.zeros(n_el)
        for el in range(0,n_el):
            No1 = int(conec_el[el,0])
            No2 = int(conec_el[el,1])
            L = math.sqrt((coord_no[No2-1,1] - coord_no[No1-1,1])**2+(coord_no[No2-1,0] - coord_no[No1-1,0])**2)
            he = L
            if (coord_no[No2-1,0] - coord_no[No1-1,0])==0:
                beta = 2*math.atan(1)
            else: beta = math.atan((coord_no[No2-1,1] - coord_no[No1-1,1])/(coord_no[No2-1,0] - coord_no[No1-1,0]))

            for i in range(0,ngl_el):
                u_el[i,0] = u[conec[el,i]-1]
            Kel = Trelica.Elem(A[el],E[el],beta,he)

            Fel[:] = np.matmul(Kel,u_el).T
            for i in range(0,2):
                Feq = math.sqrt(Fel[i]**2+Fel[i+1]**2)

            Stress[el] = Feq/A[el]

            if Fel[2]*(coord_no[No2-1,0] - coord_no[No1-1,0])<0: #F[3] corresponde a força em x no segundo nó.
                Stress[el] = -Stress[el]
        return u, Stress
