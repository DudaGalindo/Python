import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ..solver import solve_one_phase1D, solve_one_phase2D
import unittest

class test_trabalho(unittest.TestCase):
    def testA(self):
        # Entrada de Dados:
        n, P1, Pn, k = 6, 300., 50., 100
        k = k*np.ones(n)

        # Solver
        P = solve_one_phase1D.pressure(n,P1,Pn,k) #pressões tomadas no centróide do elemento

        # Validação
        def P_ans(x): return 300-250/5*x
        for x in range(n-1):
            self.assertEqual(np.round(P[x],2),P_ans(x),'ValueError')
        #plt.figure(0)
        #plt.plot(P)
        #plt.xlabel('x')
        #plt.ylabel('Pressure')
        #plt.show()

    def testB(self):
        n, P1, Pn, k1, k2 = 6,300,50,100,50
        k = np.array([k1,k1,k1,k2,k2,k2])
        P = solve_one_phase1D.pressure(n,P1,Pn,k)

        '''print(P)
        plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()'''

    def testC(self):
        n, P1, Pn, k1, k2, k3 = 6,300,50,100,50,10
        k = np.array([k1,k1,k2,k2,k3,k3])
        P = solve_one_phase1D.pressure(n,P1,Pn,k)

        #print(P)
        '''plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()'''

    def testD(self):
        #FAZER

    def testE(self):
        nx,ny, P1, Pn, k = 4, 4, 300., 50., 100
        kx = k*np.ones(nx*nx+1)
        ky = k*np.ones((ny+1)*ny)

        for i in range(ny):
            ky[i] = 0;
            ky[4*ny+i] = 0;

        for i in range(nx+1):
            kx[4*i] = 0;
        P = solve_one_phase2D.pressure(nx,ny,P1,Pn,kx,ky)
        #print(P)

        '''plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()'''

    def testeF(self):
        nx,ny,P1,Pn,k1,k2 = 4,4,300.,50.,100,1
        kx = k1*np.ones(nx*nx+1)
        ky = k1*np.ones((ny+1)*ny)

        ky[9] = k2; ky[10] = k2; ky[13] = k2; ky[14] = k2;
        ky = solve_one_phase2D.permeability(ky)

        for i in range(ny):
            ky[i] = 0;
            ky[4*ny+i] = 0;

        kx[6] = k2; kx[7] = k2; kx[10] = k2; kx[11] = k2
        kx = solve_one_phase1D.permeability(kx)

        for i in range(0,nx+1):
            kx[4*i] = 0;

        P = solve_one_phase2D.pressure(nx,ny,P1,Pn,kx,ky)
        print(P)
