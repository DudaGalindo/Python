import numpy as np
import math
import sympy as sp
from .FEM import general


class Viga:

    def init(n_el):
        n_nos_el = 2*np.ones(n_el) #número de nós por elemento
        n_nos_el = n_nos_el.astype(int)
        n_nos_tot = int(sum(n_nos_el) - (n_el-1))   #número total de nós da estrutura
        ngl_no = 2*np.ones(n_nos_tot) #numero de gl por nó
        return n_nos_el,n_nos_tot,ngl_no

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
        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no[0],n_nos_el[0],n_nos_tot)

        Fc = np.zeros(ngl_tot)
        for i in range(len(xF)):
            Fc[xF[i]] = F[i]

        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no[0])
        Kg,Fg = general.initialize(ngl_tot)
        Kg,Fg = Viga.Global(E,I,cf,Fc,x,Kg,Fg,n_el,ngl_el,conec,q)
        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg = np.linalg.inv(Kg)
        u = np.matmul(Kg,Fg.T)
        return u

    def Plots(x,u,L,n_el,E,I):

        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no[0],n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no[0])

        xe = sp.Symbol('xe'); xe1 = sp.Symbol('xe'); xe2 = sp.Symbol('xe2')

        def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
        def N2(xe,xe1,xe2):return -1*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
        def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
        def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))

        def d1N1(xe,xe1,xe2):return sp.Derivative(N1(xe,xe1,xe2),xe).doit()
        def d1N2(xe,xe1,xe2):return sp.Derivative(N2(xe,xe1,xe2),xe).doit()
        def d1N3(xe,xe1,xe2):return sp.Derivative(N3(xe,xe1,xe2),xe).doit()
        def d1N4(xe,xe1,xe2):return sp.Derivative(N4(xe,xe1,xe2),xe).doit()

        def d2N1(xe,xe1,xe2):return sp.Derivative(d1N1(xe,xe1,xe2),xe).doit()
        def d2N2(xe,xe1,xe2):return sp.Derivative(d1N2(xe,xe1,xe2),xe).doit()
        def d2N3(xe,xe1,xe2):return sp.Derivative(d1N3(xe,xe1,xe2),xe).doit()
        def d2N4(xe,xe1,xe2):return sp.Derivative(d1N4(xe,xe1,xe2),xe).doit()

        def d3N1(xe,xe1,xe2):return sp.Derivative(d2N1(xe,xe1,xe2),xe).doit()
        def d3N2(xe,xe1,xe2):return sp.Derivative(d2N2(xe,xe1,xe2),xe).doit()
        def d3N3(xe,xe1,xe2):return sp.Derivative(d2N3(xe,xe1,xe2),xe).doit()
        def d3N4(xe,xe1,xe2):return sp.Derivative(d2N4(xe,xe1,xe2),xe).doit()

        #we = np.zeros(10*n_el)
        end = 0
        n = 50
        we = np.zeros(n*n_el); thetae = np.zeros(n*n_el);Me = np.zeros(n*n_el); Ve = np.zeros(n*n_el)

        for el in range(0,n_el):
            xe1 = x[el]
            xe2 = x[el+1]
            xen = xe1

            for i in range(end,n*(el+1)):
                we[i] = -(N1(xen,xe1,xe2)*u[conec[el,0]-1]+N2(xen,xe1,xe2)*u[conec[el,1]-1]+N3(xen,xe1,xe2)*u[conec[el,2]-1]+N4(xen,xe1,xe2)*u[conec[el,3]-1])
                thetae[i] = -(d1N1(xe,xe1,xe2).subs(xe,xen)*u[conec[el,0]-1]+d1N2(xe,xe1,xe2).subs(xe,xen)*u[conec[el,1]-1]+
                d1N3(xe,xe1,xe2).subs(xe,xen)*u[conec[el,2]-1]+d1N4(xe,xe1,xe2).subs(xe,xen)*u[conec[el,3]-1])

                Me[i] = -E*I*(d2N1(xe,xe1,xe2).subs(xe,xen)*u[conec[el,0]-1] +
                d2N2(xe,xe1,xe2).subs(xe,xen)*u[conec[el,1]-1] + d2N3(xe,xe1,xe2).subs(xe,xen)*u[conec[el,2]-1] +
                d2N4(xe,xe1,xe2).subs(xe,xen)*u[conec[el,3]-1])

                Ve[i] = -E*I*(d3N1(xe,xe1,xe2).subs(xe,xen)*u[conec[el,0]-1] +
                d3N2(xe,xe1,xe2).subs(xe,xen)*u[conec[el,1]-1] + d3N3(xe,xe1,xe2).subs(xe,xen)*u[conec[el,2]-1] +
                d3N4(xe,xe1,xe2).subs(xe,xen)*u[conec[el,3]-1])

                xen = xen + (xe2-xe1)/(n-1)
            end = n*(el+1)

        ##Ploting:
        xe = np.linspace(0,L,n*n_el)

        '''plt.figure(0)
        plt.plot(xe,we)
        plt.xlabel('x')
        plt.ylabel('w')
        plt.title('Deslocamento')
        plt.show()

        plt.figure(1)
        plt.plot(xe,thetae)
        plt.xlabel('x')
        plt.ylabel('theta')
        plt.title('Rotação')
        plt.show()

        plt.figure(2)
        plt.plot(xe,Me)
        plt.xlabel('x')
        plt.ylabel('M')
        plt.title('Momento Fletor')
        plt.show()'''
        #print(Me)

        '''plt.figure(3)
        plt.plot(xe,Ve)
        plt.xlabel('x')
        plt.ylabel('V')
        plt.title('Esforço Cortante')
        plt.show()'''
        #print(Ve)
