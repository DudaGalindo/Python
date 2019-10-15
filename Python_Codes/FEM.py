import numpy as np
#from scipy.integrate import quad

class Elem:

    def Barra(n_nos_el,integralAdx,E,he,f):

        if n_nos_el == 3:
            k = 1./3*np.array([[7., -8., 1.],[-8., 16., -8.],[1., -8., 7.]])
            fel = np.array([1./6, 2./3, 1./6])

        if n_nos_el == 2:
            k = np.array([[1., -1.],[-1., 1.]])
            fel = np.array([1./2, 1./2])

        Kel = (E/he**2)*integralAdx*k    #quad(A,x1,x2)
        Fel = f*he*fel

        return Kel,Fel

    def Viga(E,I,he,cf,q,n_el):
        k1 = np.array([[6.,-3.*he,-6.,-3.*he],[-3.*he,2.*he**2,3.*he,he**2],[-6.,3.*he,6.,3*he],[-3*he,he**2,3*he,2*he**2]])
        k2 = np.array([[156.,-22*he,54.,13*he],[-22*he,4*he**2,-13*he,-3*he**2],[54.,-13*he,156.,22*he],[13*he,-3*he**2,22*he,4*he**2]])
        fel = np.array([6.,-he,6.,he])

        Kel = 2*E*I/he**3*k1 + cf*he/420*k2

        for i in range(n_el):
            q1 = q(x1,q0,L)
            q2 = q(x2,q0,L)
            f = (q1+q2)/2

        Fel = f*he/12*fel

        return Kel, Fel

class Global:
    def initialize(n_nos_tot):
        Kg = np.zeros([n_nos_tot,n_nos_tot])
        Fg = np.zeros(n_nos_tot)
        return Kg, Fg

    def Barra(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec,Fc,f,Kg,Fg):

        ##conec = np.array([[1,2],[2,3]])
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Elem.Barra(n_nos_el[el],integralAdx(he,el),E[el],he,f[el])
            for i in range(0,n_nos_el[el]):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                for j in range(0,n_nos_el[el]):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]

        Fg = Fg + Fc
        return Kg, Fg

    def Viga(E,I,cf,Fc,x,Kg,Fg,n_el,n_nos_el,conec,q):

        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Elem.Viga(E,I,he,cf,q,n_el)
            for i in range(0,n_nos_el[el]):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                for j in range(0,n_nos_el[el]):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]

        Fg = Fg + Fc
        return Kg, Fg

    def Kg_Fg_contorno(n_nos_tot,xCC,valor_CC,Kg,Fg):
        CC = np.zeros(n_nos_tot)

        for i in range(0,len(xCC)):
            CC[xCC[i]] = 1
            Fg[xCC[i]] = valor_CC[i]
            Kg[xCC[i],xCC[i]] = 1

        for i in range(0,n_nos_tot):
            if CC[i] == 1:
                for j in range(0,n_nos_tot):
                    if j!=i:
                        Kg[i,j] = 0
                        if CC[j]==0:
                            Fg[j] = Fg[j] - Kg[j,i]*Fg[i]
                            Kg[j,i] = 0

        return Kg,Fg

class deslocamento:
    def Barra(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC):
        Kg,Fg = Global.initialize(n_nos_tot)
        Kg,Fg = Global.Barra(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec,Fc,f,Kg,Fg)
        Kg,Fg = Global.Kg_Fg_contorno(n_nos_tot,xCC,valor_CC,Kg,Fg)
        Kg_inv = np.linalg.inv(Kg)
        Fg = (Fg.T)
        u = np.matmul(Kg_inv,Fg)
        return u

    def Viga(ngl_tot,n_el,n_nos_el,x,I,E,conec,Fc,q,xCC,valor_CC,cf,ngl_no):
        for i in range(0,ngl_tot):
            conec[i,0] = i
            conec[i,1] = i+1
        Kg,Fg = Global.initialize(n_nos_tot)
        Kg,Fg = Global.Viga(E,I,cf,Fc,x,Kg,Fg,n_el,n_nos_el,conec,q)
        Kg,Fg = Global.Kg_Fg_contorno(n_nos_tot,xCC,valor_CC,Kg,Fg)
        Kg_inv = np.linalg.inv(Kg)
#        Fg = (Fg.T)
#        u = np.matmul(Kg_inv,Fg)
#        return u



#class FatorForma: para vigas
#    def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
#    def N2(xe,xe1,xe2):return -3*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
#    def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
#    def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))
