import numpy as np
import matplotlib.tri as mtri
from .FEM import general


class HeatTransfer:

    def Elem_triang(K,xel,yel,fel):
        a11 = K
        a22 = K
        a12 = 0; a21 = 0; aoo = 0

        beta = np.zeros(3); alpha = np.zeros(3); gamma = np.zeros(3) #initializing
        beta[0] = yel[1] - yel[2]
        beta[1] = yel[2] - yel[0]
        beta[2] = yel[0] - yel[1]

        alpha[0] = xel[1]*yel[2] - xel[2]*yel[1]
        alpha[1] = xel[2]*yel[0] - xel[0]*yel[2]
        alpha[2] = xel[0]*yel[1] - xel[1]*yel[0]

        gamma[0] = -(xel[1] - xel[2])
        gamma[1] = -(xel[2] - xel[0])
        gamma[2] = -(xel[0] - xel[1])

        Ae = sum(alpha)/2

        beta_reshape = beta*np.ones((len(beta),len(beta)))
        gamma_reshape = gamma*np.ones((len(gamma),len(gamma)))
        Kel = (a11*beta[:,np.newaxis]*beta_reshape + a22*gamma[:,np.newaxis]*gamma_reshape)/(4*Ae)
        for i in range(0,3):
            for j in range(0,3):
                Kel[i,j] = (a11*beta[i]*beta[j]+a22*gamma[i]*gamma[j])/(4*Ae)

        Fel = fel*Ae/3*np.ones(3)
        return Kel, Fel

    def Global_triang(n_el,ngl_el,K,coord_x,coord_y,fel,conec_el_nos,Kg,Fg):

        for el in range(0,int(2*n_el)):
            No1 = int(conec_el_nos[el,0])
            No2 = int(conec_el_nos[el,1])
            No3 = int(conec_el_nos[el,2])

            xel = np.array([coord_x[No1-1],coord_x[No2-1],coord_x[No3-1]])
            yel = np.array([coord_y[No1-1],coord_y[No2-1],coord_y[No3-1]])
            Kel, Fel = HeatTransfer.Elem_triang(K,xel,yel,fel)

            for k in range(0,ngl_el):
                ig = conec_el_nos[el,k] - 1
                Fg[ig] = Fg[ig] + Fel[k]
                for p in range(0,ngl_el):
                    jg = conec_el_nos[el,p] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[k,p]

        return Kg, Fg

    def conect_el_triangular(n_el,coord_x):
        conec_el_nos = np.zeros((int(2*n_el),3))
        d=1;i=0

        for el in range(0,n_el):
            if el == (len(coord_x)-1)*d:
                d = d+1
                i = 1+i

            conec_el_nos[el*2,0] = el+1+i
            conec_el_nos[el*2,1] = el+2+i
            conec_el_nos[el*2,2] = el+len(coord_x)+1+i #ordem do mapeamento - sentido anti-horario

            conec_el_nos[el*2+1,0] = el+2+i
            conec_el_nos[el*2+1,1] = el+len(coord_x)+2+i
            conec_el_nos[el*2+1,2] = el+len(coord_x)+1+i #ordem do mapeamento - sentido anti-horario

        conec_el_nos = conec_el_nos.astype(int)
        return conec_el_nos

    def Elem_retang(K,lx,ly,fel):
        a11 = K
        a22 = K
        a12 = 0; a21 = 0; aoo = 0

        a = lx
        b = ly

        S11 = b/(6*a)*np.array([[2,-2,-1,1],[-2,2,1,-1],[-1,1,2,-2],[1,-1,-2,2]])
        S12 = 1/4*np.array([[1,1,-1,-1],[-1,-1,1,1],[-1,-1,1,1],[1,1,-1,-1]])
        S22 = a/(6*b)*np.array([[2,1,-1,-2],[1,2,-2,-1],[-1,-2,2,1],[-2,-1,1,2]])
        Soo = a*b/36*np.array([[4,2,1,2],[2,4,2,1],[1,2,4,2],[2,1,2,4]])
        Fel = 1/4*fel*a*b*np.ones(4)
        Kel = aoo*Soo + a11*S11 + a12*S12 + a22*S22

        return Kel,Fel

    def Global_retang(n_el,ngl_el,K,coord_x,coord_y,fel,conec_el_nos,Kg,Fg):
        n_nos_el = 2*np.ones(n_el).astype(int)

        for el in range(0,n_el):
            No1 = int(conec_el_nos[el,0])
            No2 = int(conec_el_nos[el,1])
            No3 = int(conec_el_nos[el,3])

            hex = abs(coord_x[No1-1] - coord_x[No2-1])
            hey = abs(coord_y[No1-1] - coord_y[No3-1])
            Kel, Fel = HeatTransfer.Elem_retang(K,hex,hey,fel)

            for k in range(0,ngl_el):
                ig = conec_el_nos[el,k] - 1
                Fg[ig] = Fg[ig] + Fel[k]
                for p in range(0,ngl_el):
                    jg = conec_el_nos[el,p] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[k,p]

        return Kg, Fg

    def conect_el(n_el,coord_x):
        conec_el_nos = np.zeros((n_el,4))
        d=1;i=0

        for el in range(0,n_el):

            if el == (len(coord_x)-1)*d:
                d = d+1
                i = 1+i

            conec_el_nos[el,0] = el+1+i
            conec_el_nos[el,1] = el+2+i
            conec_el_nos[el,2] = el+len(coord_x)+2+i #ordem do mapeamento - sentido anti-horario
            conec_el_nos[el,3] = el+len(coord_x)+1+i
        conec_el_nos = conec_el_nos.astype(int)
        return conec_el_nos

    def Temperature(n_el,x,y,coords_x,coords_y,f,xCC,valor_CC,type_element,K):
        if type_element == 'retangular':
            n_nos_el = 4*np.ones(n_el)
        else: n_nos_el = 3*np.ones(n_el)

        n_nos_tot = (len(coords_x))*(len(coords_y))
        ngl_no = 1
        ngl_tot = n_nos_tot
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)

        Kg,Fg = general.initialize(ngl_tot)

        coord_y = y[0,:]
        coord_x = coords_x
        for i in range(1,len(y[:,0])):
            coord_x = np.concatenate((coord_x,x[i,:]),axis=None)
            coord_y = np.concatenate((coord_y,y[i,:]),axis=None)

        if type_element == 'retangular':
            conec_el_nos = HeatTransfer.conect_el(n_el,coords_x)
            Kg,Fg = HeatTransfer.Global_retang(n_el,ngl_el,K,coord_x,coord_y,f,conec_el_nos,Kg,Fg)
        elif type_element == 'triangular':
            conec_el_nos = HeatTransfer.conect_el_triangular(n_el,coords_x)
            Kg,Fg = HeatTransfer.Global_triang(n_el,ngl_el,K,coord_x,coord_y,f,conec_el_nos,Kg,Fg)
        else:
            print("Not a valid entry")
            return 0

        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg_inv = np.linalg.inv(Kg)
        Fg = (Fg.T)
        T = np.matmul(Kg_inv,Fg)
        '''for el in range(n_el):
            for i in range(init,end):
                T_ans[el] = sum(T[i]*PSI(xe,ye,i))'''
        return T
