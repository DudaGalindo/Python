import numpy as np
import math
from ..FEM import general
from ..FEM_viga import Viga
from ..FEM_trelica import Trelica
import unittest

class TestTrelicas(unittest.TestCase):
    def test5(self):
        n_el = 21
        E = 200e9*np.ones(n_el)
        A = 0.0025*np.ones(n_el)
        n_nos_el = 2*np.ones(n_el)
        n_nos_tot = 12

        # Condições de Contorno:
        xCC = np.array([0,1,22]) #dos nós
        valor_CC = np.array([0.,0.,0.])

        # Forças Concentradas:
        Fc = np.array([-20E3,-30E3,-40E3])
        xFc = np.array([2,6,10])

        coord_no  = np.array([[0.,0.],[2.,0.],[2.,3.],[4.,0.],[4.,3.],[6.,0.],[6.,3.],\
                              [8.,3.],[8.,0.],[10.,3.],[10.,0],[12.,0.]]) #[x,y]

        conec_el = np.array([[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6],\
                             [5,7],[6,7],[6,8],[6,9],[7,8],[8,9],[8,10],[9,10],[9,11],\
                             [10,11],[10,12],[11,12]])

        u,T = Trelica.deslocamento_tensao(coord_no,conec_el,n_nos_tot,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E)

        #print(u)
        #print(T)

    def test6(self):
        n_el = 17
        E = 200e9*np.ones(n_el)
        A = 0.0025*np.ones(n_el)
        n_nos_el = 2*np.ones(n_el)
        n_nos_tot = 10

        # Condições de Contorno:
        xCC = np.array([0,18]) #dos nós
        valor_CC = np.array([0.,0.])

        # Forças Concentradas:
        Fc = np.array([80,80,60,40])*4.448  #[N]
        xFc = np.array([2,6,10,14])

        coord_no  = np.array([[0.,0.],[2.,0.],[2.,2.],[4.,0.],[4.,2.],[6.,0.],[6.,2.],\
                              [7.5,0],[7.5,2],[0,2.]])*0.305 #[x,y] [m]
        conec_el = np.array([[1,2],[1,3],[1,10],[10,3],[2,3],[2,4],[2,5],[3,5],[4,5],\
                            [4,6],[4,7],[5,7],[6,7],[6,8],[6,9],[7,9],[8,9]])
        u,T = Trelica.deslocamento_tensao(coord_no,conec_el,n_nos_tot,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E)

class TestVigas(unittest.TestCase):
    def test3_general_case(self): #checar os dados de entrada
        n_el = 50 #numero de elementos da estrutura
        L = 13*0.305 #m
        E = 200e9 #constante
        I = 50e-6 #constante

        def q(x):
            q0 = 250*4.448/0.305
            q1 = 150*4.448/0.305
            if x<8*0.305:
                return q0*x/(8*0.305)
            else: return q0+q1*(x-8*0.305)/(5*0.305)

        # Forças Concentradas:
        n_el_trechos = np.array([round(8/13*n_el),round(2/13*n_el),round(3/13*n_el)])
        xFc = np.array([2*n_el_trechos[0],2*(n_el_trechos[0]+n_el_trechos[1])])
        Fc = np.array([3000,-15000])*4.448

        cf = 0;

        # Condições de Contorno:
        n_nos_el,n_nos_tot,ngl_no = Viga.init(n_el)
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el[0],n_nos_tot)

        xCC_dir = np.array([ngl_tot-2,ngl_tot-1]) #posição
        valor_CCdir = np.array([0, 0])

        # Vetor x:
        x1 = np.linspace(0,8*0.305,n_el_trechos[0]+1)
        x2 = np.linspace(8.2*0.305,10*0.305,n_el_trechos[1])
        x3 = np.linspace((10.2)*0.305,L,n_el_trechos[2])
        x = np.append(x1,x2); x = np.append(x,x3)

        u = Viga.deslocamento(n_el,x,I,E,Fc,xFc,q,xCC_dir,valor_CCdir,cf) #1GL
        Viga.Plots(x,u,L,n_el,E,I)

class testeFRAME(unittest.TestCase):
    def test
