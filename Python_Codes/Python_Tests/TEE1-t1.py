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
        xPn = np.array([1,n])
        # Solver
        P = solve_one_phase1D.pressure(n,P1,Pn,xPn,k) #pressões tomadas no centróide do elemento

        # Validação
        def P_ans(x): return 300-250/5*x
        for x in range(n-1):
            self.assertEqual(np.round(P[x],2),P_ans(x),'ValueError')
        print('A:',P)
        #plt.figure(0)
        #plt.plot(P)
        #plt.xlabel('x')
        #plt.ylabel('Pressure')
        #plt.show()

    def testB(self):
        n, P1, Pn, k1, k2 = 6,300,50,100,50
        k = np.array([k1,k1,k1,k2,k2,k2])
        xPn = np.array([1,n])
        P = solve_one_phase1D.pressure(n,P1,Pn,xPn,k)
        print('B:',P)

        plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()

    def testC(self):
        n, P1, Pn, k1, k2, k3 = 6,300,50,100,50,10
        k = np.array([k1,k1,k2,k2,k3,k3])
        xPn = np.array([1,n])
        P = solve_one_phase1D.pressure(n,P1,Pn,xPn,k)
        print('C:',P)

        plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()

    def testD(self):
        nx,P1,qn,k,L = 6, 300.,250.,100,100
        qNewman = 250
        xNewman = nx
        xP = np.array([1])
        k = k*np.ones(nx);
        P = solve_one_phase1D.pressure_Newmann(nx,P1,xP,qNewman,xNewman,k,L)
        print('D:',P)

    def testE(self):
        nx,ny, P1, Pn, k = 4, 4, 300., 50., 100
        k = k*np.ones([nx,ny])

        P = solve_one_phase2D.pressure(nx,ny,P1,Pn,k)
        print('E:',P)

    def testeF(self):
        nx,ny,P1,Pn,k1,k2 = 4,4,300.,50.,100,1
        k = k1*np.ones([nx,ny])
        k[1,1] = k2; k[2,2] = k2; k[1,2] = k2; k[2,1] = k2;

        P = solve_one_phase2D.pressure(nx,ny,P1,Pn,k)
        print('F:',P)
