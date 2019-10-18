import numpy as np
import math
from ..FEM import Viga, Barra, general
import matplotlib.pyplot as plt
import unittest

class Test_FEM(unittest.TestCase):
    def testFEM_Barra_prob5_4(self): # 1grau de liberdade por nó - cálculo do deslocamento apenas
    # ENTRADA DE DADOS
        n_el = 2 #numero de elementos da estrutura (AB e BC)
        n_nos_el = np.array([3,2]) #número de nós por elemento
        ngl_no = 1
        x = np.array([1., 3., 5.]); xA = x[0]; xB = x[1]; xC = x[2]
        n_nos_tot = int(sum(n_nos_el) - (n_el-1)) #número total de nós da estrutura

        E = 2E7*np.ones(n_el)

        Fc = np.array([0, 0, 150., 0])
        f = np.array([10.,0])

        def integralAdx(he,i):
            if i==0: return 0.1*he
            else: return 0.3

        conec = np.array([[1,2,3], [3,4,0]])    ##barra AB tem 3nós e a BC tem 2. - caso tenha alguma outra coisa doida assim, preferi deixar logo pra o usuario definir

        # posição e valor das condições de contorno
        xCC = np.array([0,n_nos_tot-1]) #posição
        valor_CC = np.array([0, 0])

        u = Barra.deslocamento(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC) #1GL

        u_ans = np.array([0.,3.45E-5,6.4E-5,0.])
        for i in range(len(u)):
            self.assertAlmostEqual(u[i],u_ans[i],3,'ValueError: Failed')

    def testFEM_Viga(self): #checar os dados de entrada se possível
        n_el = 8 #numero de elementos da estrutura
        L = 2.
        E = 200e9 #constante
        I = 50e-6 #constante
        q0 = 200E3 #cte

        def q(x):
            q0 = 200E3
            L = 2.
            return q0*np.sin((2*math.pi*x)/L)

        x = np.linspace(0,2.,n_el+1);

        Fc = np.array([-50E3,50E3])
        xFc = np.array([4,12])

        cf = 0;

        xCC_dir = np.array([0,1,16,17]) #posição
        valor_CCdir = np.array([0, 0, 0, 0])
        u = Viga.deslocamento(n_el,x,I,E,Fc,xFc,q,xCC_dir,valor_CCdir,cf) #1GL
        print(u)
        u_ans = np.array([0.,0.,1.5025E-5,-8.4188E-5,3.0318E-5,3.3573E-5,2.7515E-5,6.7402E-5,1.7364E-20,1.3429E-4,-2.7615E-5,6.7402E-5,-3.0318E-5,-3.3573E-5,-1.5025E-5,-8.4188E-5,0.,0.])
        for i in range(len(u)):
            self.assertAlmostEqual(u[i],u_ans[i],3,'ValueError: Failed')

    #    xe = np.linspace(0,2.,10*(n_el  ))
        we = Viga.Plots(x,u,L,n_el)
        print(we)
        xe = np.linspace(0,L,len(we))
        plt.figure(0)
        plt.plot(xe,we)
        plt.xlabel('x')
        plt.ylabel('we')
        plt.show()
