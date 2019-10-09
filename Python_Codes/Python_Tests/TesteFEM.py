import numpy as np
from ..FEM import Kelem
import unittest

class Test_FEM(unittest.TestCase):
    def FEM_1GL(self): # 1grau de liberdade por nó - cálculo do deslocamento apenas
    # ENTRADA DE DADOS
        n_el = 2; #numero de elementos da estrutura
        n_nos_AB = 3; #número de nós em cada elemento
        n_nos_BC = 2;
        n_nos_tot = n_nos_AB*1 + n_nos_BC*1-(n_el-1); #número total de nós da estrutura
        A = np.array([0.1, 0.1, 100])
        E = np.array([200, 200, 200])
        he = np.array([])
        Ngl_no = 1;

    # INICIALIZAÇÃO:
        for i in range(n_el):
            Kel = Kelem.Kel(n_nos_el,A(i),E(i),he(i))
