import numpy as np
import math
from ..FEM import general
from ..FEM_frame import Frame
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
        n_el_trechos = (np.array([8/13,2/13,3/13])*n_el).round().astype(int)
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
        x2 = np.linspace((8+2/n_el_trechos[1])*0.305,10*0.305,n_el_trechos[1])
        x3 = np.linspace((10+3/n_el_trechos[2])*0.305,L,n_el_trechos[2])
        x = np.append(x1,x2); x = np.append(x,x3)

        u = Viga.deslocamento(n_el,x,I,E,Fc,xFc,q,xCC_dir,valor_CCdir,cf) #1GL
        Viga.Plots(x,u,L,n_el,E,I)

class testeFRAME(unittest.TestCase):
    def test6(self):
        E = 30E6 #lb/in²
        n_el = 104
        A = 100*np.ones(n_el) #in²
        I = 200*np.ones(n_el) #in⁴

        n_nos_el = 2*np.ones(n_el)
        coord_no = np.array([[0,0],[0,16],[8,16],[20,16],[20,0]])*12

        def q(x):
            if x<16*12:return -500
            else: return 0

        n_el_trechos = (np.array([16/52,8/52,12/52,16/52])*n_el).round().astype(int)

        coord_no_trecho1 = np.linspace(coord_no[0,:],coord_no[1,:],n_el_trechos[0]+1)
        coord_no_trecho2 = np.linspace(coord_no[1,:]+np.array([8*12/n_el_trechos[1],0]),coord_no[2,:],n_el_trechos[1])
        coord_no_trecho3 = np.linspace(coord_no[2,:]+np.array([12*12/n_el_trechos[2],0]),coord_no[3,:],n_el_trechos[2])
        coord_no_trecho4 = np.linspace(coord_no[3,:]-np.array([0,16*12/n_el_trechos[3]]),coord_no[4,:],n_el_trechos[3])
        coord_nox = np.append(coord_no_trecho1[:,0],coord_no_trecho2[:,0])
        coord_nox = np.append(coord_nox[:],coord_no_trecho3[:,0])
        coord_nox = np.append(coord_nox[:],coord_no_trecho4[:,0])
        coord_noy = np.append(coord_no_trecho1[:,1],coord_no_trecho2[:,1])
        coord_noy = np.append(coord_noy[:],coord_no_trecho3[:,1])
        coord_noy = np.append(coord_noy[:],coord_no_trecho4[:,1])
        coord_no = np.array([coord_nox,coord_noy]).T

        # Esforços concentrados/prescritos:
        Fc = np.array([-10E3]) #checar isso e o gl de aplicação
        xFc = np.array([3*(n_el_trechos[0]+n_el_trechos[1])-1-1])

        # Condições de Contorno:
        n_nos_tot = general.n_nosTOTAL(n_nos_el,n_el)
        ngl_no = 3
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el,n_nos_tot)

        xCC = np.array([0,1,ngl_tot-3,ngl_tot-2,ngl_tot-1])
        valor_CC = np.array([0,0,0,0,0])
        u,reacoes = Frame.deslocamento_reacao(coord_no,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E,I,q)


    def test7(self):
        E = 30E6 #lb/in²
        n_el = 120
        A = 100*np.ones(n_el) #in²
        I = 200*np.ones(n_el) #in⁴

        n_nos_el = 2*np.ones(n_el)
        coord_no = np.array([[0,0],[0,10],[10,10],[10,0]])*12 #12 é conversão de ft pra in

        def q(x):
            return 0

        n_el_trechos = (np.array([10/30,10/30,10/30])*n_el).round().astype(int)

        for i in range(0,n_el_trechos[0]):
            I[i+n_el_trechos[0]] = 100


        coord_no_trecho1 = np.linspace(coord_no[0,:],coord_no[1,:],n_el_trechos[0]+1)
        coord_no_trecho2 = np.linspace(coord_no[1,:]+np.array([10*12/n_el_trechos[1],0]),coord_no[2,:],n_el_trechos[1])
        coord_no_trecho3 = np.linspace(coord_no[2,:]-np.array([0,10*12/n_el_trechos[2]]),coord_no[3,:],n_el_trechos[2])
        coord_nox = np.append(coord_no_trecho1[:,0],coord_no_trecho2[:,0])
        coord_nox = np.append(coord_nox[:],coord_no_trecho3[:,0])
        coord_noy = np.append(coord_no_trecho1[:,1],coord_no_trecho2[:,1])
        coord_noy = np.append(coord_noy[:],coord_no_trecho3[:,1])
        coord_no = np.array([coord_nox,coord_noy]).T


        # Esforços concentrados/prescritos:
        Fc = np.array([10E3,-5000]) #checar isso e o gl de aplicação
        xFc = np.array([3*n_el_trechos[0]-1+1,3*(n_el_trechos[0]+n_el_trechos[1])-1])

        # Condições de Contorno:
        n_nos_tot = general.n_nosTOTAL(n_nos_el,n_el)
        ngl_no = 3
        ngl_el, ngl_tot = general.ngl(ngl_no,n_nos_el,n_nos_tot)

        xCC = np.array([0,1,ngl_tot-3,ngl_tot-2,ngl_tot-1])
        valor_CC = np.array([0,0,0,0,0])
        u,reacoes = Frame.deslocamento_reacao(coord_no,n_nos_el,Fc,xFc,xCC,valor_CC,n_el,A,E,I,q)
