import numpy as np
#from scipy.integrate import quad

class Elem:

    def Kel_Fel(n_nos_el,integralAdx,E,he,f):
        if n_nos_el == 3:
            k = 1./3*np.array([[7., -8., 1.],[-8., 16., -8.],[1., -8., 7.]])
            fel = np.array([1./6, 2./3, 1./6])

        if n_nos_el == 2:
            k = np.array([[1., -1.],[-1., 1.]])
            fel = np.array([1./2, 1./2])

        Kel = (E/he**2)*integralAdx*k    #quad(A,x1,x2)
        Fel = f*he*fel

        return Kel,Fel

class Global:
    def Kg_Fg(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec,Fc,f):
        Kg = np.zeros([n_nos_tot,n_nos_tot])
        Fg = np.zeros(n_nos_tot)

        ##conec = np.array([[1,2],[2,3]])
        for el in range(0,n_el):
            he = abs(x[el] - x[el+1])
            Kel, Fel = Elem.Kel_Fel(n_nos_el[el],integralAdx(he,el),E[el],he,f[el])

            for i in range(0,n_nos_el[el]):
                ig = conec[el,i] - 1
                Fg[ig] = Fg[ig] + Fel[i]
                #if el>0 and i == 0: ig = ig
                for j in range(0,n_nos_el[el]):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]

        Fg = Fg + Fc
        return Kg, Fg


    def Kg_Fg_contorno(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC):
        CC = np.zeros(n_nos_tot)
        Kg, Fg = Global.Kg_Fg(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec,Fc,f)

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

class ans:
    def deslocamento(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC):
        Kg,Fg = Global.Kg_Fg_contorno(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC)
        Kg_inv = np.linalg.inv(Kg)
        Fg = (Fg.T)
        u = np.matmul(Kg_inv,Fg)
        return u



#class FatorForma: para vigas
#    def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
#    def N2(xe,xe1,xe2):return -3*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
#    def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
#    def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))
