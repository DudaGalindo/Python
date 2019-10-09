import numpy as np
#from scipy.integrate import quad

class Elem:
    def hel(x1,x2,n_nos_el):
        hel = (x2-x1)/n_nos_el
        return hel

    def Kel(n_nos_el,integralAdx,E,he):
        if n_nos_el == 3:
            k = 1./3*np.array([[7., -8., 1.],[-8., 16., -8.],[1., -8., 7.]])
        if n_nos_el == 2:
            k = np.array([[1., -1.],[-1., 1.]])
        Kel = (E/he**2)*integralAdx*k    #quad(A,x1,x2)
        return Kel

    def Kg(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec):
        Kg = np.zeros([n_nos_tot,n_nos_tot])
        ##conec = np.array([[1,2],[2,3]])
        for el in range(0,n_el):
            he = Elem.hel(x[el], x[el+1], n_nos_el[el]-1)
            Kel = Elem.Kel(n_nos_el[el],integralAdx(he,el),E[el],he)
            for i in range(0,n_nos_el[el]):
                ig = conec[el,i] - 1
                #if el>0 and i == 0: ig = ig
                for j in range(0,n_nos_el[el]):
                    jg = conec[el,j] - 1
                    Kg[ig,jg] = Kg[ig,jg] + Kel[i,j]
        return Kg

    def Fel(n_nos_tot):
        for i in range(n_nos_el):
            if n_nos_el == 3: aa
            if n_nos_el == 2:
                Fel = f*he/2*[1, 1]


#class FatorForma:
#    def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
#    def N2(xe,xe1,xe2):return -3*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
#    def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
#    def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))
