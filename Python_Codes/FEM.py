import numpy as np
import math
#from scipy.integrate import quad

class general:

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

    def conect(ngl_tot,n_el,ngl_el):
        vetor_ngl = np.linspace(1,ngl_tot,ngl_tot)
        conec = np.zeros((n_el,ngl_el))
        for el in range(0,n_el):
            for i in range(0,ngl_el):
                conec[el,i] = vetor_ngl[i+el*2]
        conec = conec.astype(int)
        return conec

    def Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg):
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

    def deslocamento(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC):

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


class Viga:


    def Elem(E,I,he,cf,q,x1,x2):
        k1 = np.array([[6.,-3.*he,-6.,-3.*he],[-3.*he,2.*he**2,3.*he,he**2],[-6.,3.*he,6.,3*he],[-3*he,he**2,3*he,2*he**2]])
        k2 = np.array([[156.,-22*he,54.,13*he],[-22*he,4*he**2,-13*he,-3*he**2],[54.,-13*he,156.,22*he],[13*he,-3*he**2,22*he,4*he**2]])
        fel = np.array([6.,-he,6.,he])

        Kel = 2*E*I/(he**3)*k1 + cf*he/420*k2

        q1 = q(x1)
        q2 = q(x2)
        f = (q1+q2)/2
        Fel = f*he/12*fel

        return Kel, Fel

    def Global(E,I,cf,Fc,x,Kg,Fg,n_el,ngl_el,conec,q):

        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Viga.Elem(E,I,he,cf,q,x[el],x[el+1])
            for i in range(0,ngl_el):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                for j in range(0,ngl_el):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]

        Fg = Fg + Fc
        return Kg, Fg

    def deslocamento(n_el,x,I,E,F,xF,q,xCC,valor_CC,cf):
        n_nos_el = 2*np.ones(n_el) #número de nós por elemento
        n_nos_el = n_nos_el.astype(int)
        n_nos_tot = int(sum(n_nos_el) - (n_el-1))   #número total de nós da estrutura
        ngl_no = 2*np.ones(n_nos_tot) #numero de gl por nó
        ngl_el, ngl_tot = general.ngl(ngl_no[0],n_nos_el[0],n_nos_tot)

        Fc = np.zeros(ngl_tot)
        Fc[xF[0]] = F[0]; Fc[xF[1]] = F[1]

        conec = general.conect(ngl_tot,n_el,ngl_el)
        Kg,Fg = general.initialize(ngl_tot)
        Kg,Fg = Viga.Global(E,I,cf,Fc,x,Kg,Fg,n_el,ngl_el,conec,q)
        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg = np.linalg.inv(Kg)
        u = np.matmul(Kg,Fg.T)
        return u

    def FatorForma(x):
        def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
        def N2(xe,xe1,xe2):return -1*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
        def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
        def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))
            
        plt.figure(0)
        plt.plot()
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()

    def Fletor():

    def Cortante():
