import numpy as np
import math
import sympy as sp
import matplotlib.pyplot as plt
from .FEM import general


class Viga:

    def init(n_el):
        n_nos_el = 2*np.ones(n_el) #número de nós por elemento
        n_nos_el = n_nos_el.astype(int)
        n_nos_tot = int(sum(n_nos_el) - (n_el-1))   #número total de nós da estrutura
        ngl_no = 2 #numero de gl por nó
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
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no)

        Fc = np.zeros(ngl_tot)
        Fc[xF[:]] = F[:]

        Kg,Fg = general.initialize(ngl_tot)
        Kg,Fg = Viga.Global(E,I,cf,Fc,x,Kg,Fg,n_el,ngl_el,conec,q)
        Kg,Fg = general.Kg_Fg(ngl_tot,xCC,valor_CC,Kg,Fg)
        Kg = np.linalg.inv(Kg)
        u = np.matmul(Kg,Fg.T)
        return u

    def Plots(x,u,L,n_el,E,I):

        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)
        conec = general.conect(ngl_tot,n_el,ngl_el,ngl_no)
        xe = sp.Symbol('xe'); xe1 = sp.Symbol('xe'); xe2 = sp.Symbol('xe2')

        def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
        def N2(xe,xe1,xe2):return -1*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
        def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
        def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))

        def d1N1(xe,xe1,xe2):return -6*(xe-xe1)/((xe2-xe1)**2)+6*(xe-xe1)**2/(xe2-xe1)**3
        def d1N2(xe,xe1,xe2):return -1*(1+3*((xe-xe1)/(xe2-xe1))**2-4*((xe-xe1)/(xe2-xe1)))
        def d1N3(xe,xe1,xe2):return (6*(xe-xe1)/(xe2-xe1)**2)*(1-(xe-xe1)/(xe2-xe1))
        def d1N4(xe,xe1,xe2):return -1*((xe-xe1)/(xe2-xe1))*((3*((xe-xe1)/(xe2-xe1)))-2)

        def d2N1(xe,xe1,xe2):return (-6/(xe2-xe1)**2)*(1-2*(xe-xe1)/(xe2-xe1))
        def d2N2(xe,xe1,xe2):return -6*((xe-xe1)/(xe2-xe1)**2)+4/(xe2-xe1)
        def d2N3(xe,xe1,xe2):return 6/(xe2-xe1)**2-12*((xe-xe1)/(xe2-xe1)**3)
        def d2N4(xe,xe1,xe2):return 2/(xe2-xe1)-6*((xe-xe1)/(xe2-xe1)**2)

        def d3N1(xe,xe1,xe2):return 12/(xe2-xe1)**3
        def d3N2(xe,xe1,xe2):return -6/(xe2-xe1)**2
        def d3N3(xe,xe1,xe2):return -12/(xe2-xe1)**3
        def d3N4(xe,xe1,xe2):return -6/(xe2-xe1)**2

        #we = np.zeros(10*n_el)
        end = 0
        n = 50
        we = np.zeros(n*n_el); thetae = np.zeros(n*n_el);Me = np.zeros(n*n_el); Ve = np.zeros(n*n_el)

        for el in range(0,n_el):
            xe1 = x[el]
            xe2 = x[el+1]
            xe = np.linspace(xe1,xe2,n)

            we[end:n*(el+1)] = -(N1(xe,xe1,xe2)*u[conec[el,0]-1]+N2(xe,xe1,xe2)*u[conec[el,1]-1]+N3(xe,xe1,xe2)*u[conec[el,2]-1]+N4(xe,xe1,xe2)*u[conec[el,3]-1])

            thetae[end:n*(el+1)] = -(d1N1(xe,xe1,xe2)*u[conec[el,0]-1]+d1N2(xe,xe1,xe2)*u[conec[el,1]-1]+
            d1N3(xe,xe1,xe2)*u[conec[el,2]-1]+d1N4(xe,xe1,xe2)*u[conec[el,3]-1])

            Me[end:n*(el+1)] = -E*I*(d2N1(xe,xe1,xe2)*u[conec[el,0]-1] +
            d2N2(xe,xe1,xe2)*u[conec[el,1]-1] + d2N3(xe,xe1,xe2)*u[conec[el,2]-1] +
            d2N4(xe,xe1,xe2)*u[conec[el,3]-1])

            Ve[end:n*(el+1)] = -E*I*(d3N1(xe,xe1,xe2)*u[conec[el,0]-1] +
            d3N2(xe,xe1,xe2)*u[conec[el,1]-1] + d3N3(xe,xe1,xe2)*u[conec[el,2]-1] +
            d3N4(xe,xe1,xe2)*u[conec[el,3]-1])

            end = n*(el+1)

        ##Ploting:
        xe = np.linspace(0,L,n*n_el)

        plt.figure(0)
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
        plt.show()
        #print(Me)

        plt.figure(3)
        plt.plot(xe,Ve)
        plt.xlabel('x')
        plt.ylabel('V')
        plt.title('Esforço Cortante')
        plt.show()
        #print(Ve)
