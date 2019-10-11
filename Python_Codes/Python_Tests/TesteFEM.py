import numpy as np
from ..FEM import Elem, Global, ans
import unittest

class Test_FEM(unittest.TestCase):
    def testFEM_1GL(self): # 1grau de liberdade por nó - cálculo do deslocamento apenas
    # ENTRADA DE DADOS
        n_el = 2 #numero de elementos da estrutura (AB e BC)
        n_nos_el = np.array([3,2]) #número de nós por elemento
        x = np.array([1., 3., 5.]); xA = x[0]; xB = x[1]; xC = x[2]
        n_nos_tot = sum(n_nos_el) - (n_el-1) #número total de nós da estrutura

        E = 2e7*np.ones(n_el)

        Fc = np.array([0, 0, 150., 0])
        f = np.array([10.,0])

        def integralAdx(he,i):
            if i==0: return 0.1*he
            else: return 0.3

        conec = np.array([[1,2,3], [3,4,0]])    ##barra AB tem 3nós e a BC tem 2. - caso tenha alguma outra coisa doida assim, preferi deixar logo pra o usuario definir

        # posição e valor das condições de contorno
        xCC = np.array([0,n_nos_tot-1]) #posição
        valor_CC = np.array([0, 0])

        u = ans.deslocamento(n_nos_tot,n_el,n_nos_el,x,integralAdx,E,conec,Fc,f,xCC,valor_CC)

        u_ans = np.array([0.,3.45E-5,6.4E-5,0.])
        for i in range(len(u)):
            self.assertAlmostEqual(u[i],u_ans[i],3,'ValueError: Failed')
