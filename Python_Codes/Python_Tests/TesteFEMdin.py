''' Testes Viga Dinâmicos '''

import numpy as np
from ..FEM_viga_dinamica import VigaDin
import unittest

class TestVigasDinamica(unittest.TestCase):
    def test6(self):
        alfa = 0.000000001
        beta = 0.000000005
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

        VigaDin.deslocamento(E,I,cf,x,xCC_dir,valor_CCdir,xFc,Fc,n_el,q,alpha,beta)
