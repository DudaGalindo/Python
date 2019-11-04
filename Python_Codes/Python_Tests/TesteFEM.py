import numpy as np
import math
from ..FEM import Viga, Barra, general, Trelica, Frame
import matplotlib.pyplot as plt
import unittest

class Test_FEM(unittest.TestCase):

    '''def testFEM_Barra_prob5_4(self): # 1grau de liberdade por nó - cálculo do deslocamento apenas
    # ENTRADA DE DADOS
        n_el = 2 #numero de elementos da estrutura (AB e BC)
        n_nos_el = np.array([3,2]) #número de nós por elemento
        ngl_no = 1
        x = np.array([1., 3., 5.]);

        E = 2E7*np.ones(n_el)

        Fc = np.array([0, 0, 150., 0])
        f = np.array([10.,0])

        def integralAdx(he,i):
            if i==0: return 0.1*he
            else: return 0.3

        conec = np.array([[1,2,3], [3,4,0]])    ##barra AB tem 3nós e a BC tem 2. - caso tenha alguma outra coisa doida assim, preferi deixar logo pra o usuario definir

        # posição e valor das condições de contorno
        xCC = np.array([0,3]) #posição
        valor_CC = np.array([0, 0])

        u = Barra.deslocamento(n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC) #1GL

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

        # Forças Concentradas
        Fc = np.array([-50E3,50E3])
        xFc = np.array([4,12])

        cf = 0;

        # Condições de Contorno:
        xCC_dir = np.array([0,1,16,17]) #posição
        valor_CCdir = np.array([0, 0, 0, 0])
        x = np.linspace(0,L,n_el+1);

        u = Viga.deslocamento(n_el,x,I,E,Fc,xFc,q,xCC_dir,valor_CCdir,cf) #1GL

        u_ans = np.array([0.,0.,1.5025E-5,-8.4188E-5,3.0318E-5,3.3573E-5,2.7515E-5,6.7402E-5,1.7364E-20,1.3429E-4,-2.7615E-5,6.7402E-5,-3.0318E-5,-3.3573E-5,-1.5025E-5,-8.4188E-5,0.,0.])

        for i in range(len(u)):
            self.assertAlmostEqual(u[i],u_ans[i],3,'ValueError: Failed')
        ##Plotagens:
        Viga.Plots(x,u,L,n_el,E,I)

    def test_trelica7_4_1(self):
        n_el = 2;
        E = 3E7*np.ones(n_el)
        A = np.array([0.4,0.5])

        n_nos_el = 2*np.ones(n_el)
        n_nos_tot = general.n_nosTOTAL(n_nos_el,n_el)

        # Condições de Contorno:
        xCC = np.array([0,1,4,5]) #dos nós
        valor_CC = np.array([0.,0.,0.,0.])

        # Forças Concentradas:
        Fc = np.array([-1E3])
        xFc = np.array([3])

        coord_no  = np.array([[0.,0.],[10.,0.],[0.,10.]]) #[x,y]
        conec_el = np.array([[1,2],[2,3]])
        u = Trelica.deslocamento(coord_no,conec_el,n_nos_tot,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E)

        u_ans = np.array([0.0000,0.0000,-0.0008,-0.0027,0.0000,0.0000])

        for i in range(len(u)):
            self.assertAlmostEqual(u[i],u_ans[i],4,'ValueError: Failed')

    def test_trelica7_4_2(self): ##rever
        n_el = 9
        E = 200e9*np.ones(n_el)
        A = 0.0025*np.ones(n_el)

        n_nos_el = 2*np.ones(n_el)
        n_nos_tot = 6

        # Forças concentradas:
        Fc = np.array([-600,200])
        xFc = np.array([7,8])

        # Condições de Contorno:
        xCC = np.array([0,1,11])
        valor_CC = np.array([0,0,0])

        coord_no = np.array([[0,0],[4,0],[4,3],[8,0],[8,3],[12,0]])
        conec_el = np.array([[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6]])

        u = Trelica.deslocamento(coord_no,conec_el,n_nos_tot,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E)

        u_ans = np.array([0.0000,0.0000E-5,0.3200E-5,-1.57E-5,0.8650E-5,-1.5700E-5,0.6400E-5,-2.2867E-5,0.5450E-5,-2.0167E-5,1.12000E-5,0.0000])

        for i in range(len(u)-1):
            self.assertAlmostEqual(u[i],u_ans[i],9,'ValueError: Failed')'''

    def testeFRAME(self):
        E = 210000
        A = 0.01
        I = 2E8
        n_el = 3
        n_nos_el = 2*np.ones(n_el)
        n_nos_el = n_nos_el.astype(int)
        coord_no = np.array([[0,3],[3,3],[6,0],[9,0]])

        # Esforços concentrados/prescritos:
        Fc = np.array([-10E3,-5E3,-10E3,5E6])
        xFc = np.array([4,5,7,8])

        # Condições de Contorno:
        xCC = np.array([0,1,2,9,10,11])
        valor_CC = np.array([0,0,0,0,0,0])
        u = Frame.deslocamento(coord_no,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E,I)
        print(u)
