import numpy as np
from ..FEM import Elem
import unittest

class Test_FEM(unittest.TestCase):
    def testFEM_1GL(self): # 1grau de liberdade por nó - cálculo do deslocamento apenas
    # ENTRADA DE DADOS
        n_el = 2; #numero de elementos da estrutura (AB e BC)
        n_nos_el = np.array([3,2]); #número de nós por elemento
        x = np.array([1., 3., 5.]); xA = x[0]; xB = x[1]; xC = x[2]
        n_nos_tot = sum(n_nos_el) - (n_el-1); #número total de nós da estrutura

        E = 2e7*np.ones(n_el);

        def integralAdx(he,i):
            if i==0: return 0.1*he
            else: return 0.3

        conec = np.array([[1,2,3], [3,4,0]])    ##barra AB tem 3nós e a BC tem 2. - caso tenha alguma outra coisa doida assim, preferi deixar logo pra o usuario definir
    # INICIALIZAÇÃO
        Kg = Elem.Kg(n_el,n_nos_tot,n_nos_el,x,integralAdx,E,conec)
        print(Kg)
