import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ..solver import solve_one_phase
import unittest

class test_trabalho(unittest.TestCase):
    def testA(self):
        # Entrada de Dados:
        n, P1, Pn, k = 6, 300., 50., 100
        k = k*np.ones(n)
        
        # Solver
        P = solve_one_phase.pressure(n,P1,Pn,k) #pressões tomadas no centróide do elemento

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
        P = solve_one_phase.pressure(n,P1,Pn,k)

        #print(P)
        #plt.figure(0)
        #plt.plot(P)
        #plt.xlabel('x')
        #plt.ylabel('Pressure')
        #plt.show()

    def testC(self):
        n, P1, Pn, k1, k2, k3 = 6,300,50,100,50,10
        k = np.array([k1,k1,k2,k2,k3,k3])
        P = solve_one_phase.pressure(n,P1,Pn,k)

        print(P)
        plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()
