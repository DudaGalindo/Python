import numpy as np
import matplotlib.pyplot  as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib
from ..FEM_HeatTransfer import HeatTransfer
from scipy.spatial import Delaunay
import unittest

class test_trabalho(unittest.TestCase):
    def test_retangular_element(self):
        ''' Dados de Entrada '''
        n_elem = 49
        nx = 7 #número de elementos em x
        lx = 1
        ly = 1
        type_element = 'triangular'
        f = 0
        K_steel = 52

        # Creating mesh:
        ny = round(n_elem/nx) #número de elementos em y
        coord_x = np.linspace(0,lx,nx+1)
        coord_y = np.linspace(0,ly,ny+1)
        x,y = np.meshgrid(coord_x,coord_y)

        ''' Condições de Contorno '''
        T_CCx = np.linspace(0,nx,nx+1)
        T_CClateralesq = np.zeros(ny)
        for i in range(0,ny):
            T_CClateralesq[i] = (i+1)*(nx+1)

        T_CC = np.concatenate((T_CCx,T_CClateralesq),axis=None).astype(int)

        valorT_CC = 100*np.ones(nx+1)
        valorT_CClateralesq = 50*np.ones(ny)
        valorT_CC = np.concatenate((valorT_CC,valorT_CClateralesq),axis=None)

        T = HeatTransfer.Temperature(n_elem,x,y,coord_x,coord_y,f,T_CC,valorT_CC,type_element,K_steel)

        Treshape = np.zeros(x.shape)
        m = 0
        for ti in range(ny+1):
            i = (nx+1)*m
            for tj in range(nx+1):
                Treshape[ti,tj] = T[i+tj]
            m = m+1
        #print('resposta:',Treshape)
        print('resposta',Treshape)
        if type_element == 'retangular':
            plt.pcolormesh(x,y,Treshape,cmap = cm.RdGy)
            plt.colorbar()
            plt.show()
        else:
            plt.pcolormesh(x,y,Treshape,cmap = cm.RdGy)
            plt.colorbar()
            plt.show()
            matplotlib.axes.Axes.tripcolor(x,y,Treshape,cmap = cm.RdGy)
