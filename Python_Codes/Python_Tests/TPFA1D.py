import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ..solver import solve_one_phase, solve_two_phases
import unittest

class test_TPFA(unittest.TestCase):
    def test1D_Dirichlet(self):
        # Entrada de Dados:
        n, P1, P5, k = 5, 1., 0.,1
        #x,y = np.meshgrid(range(n),range(1))
        #plt.imshow(x)
        k = 1*np.ones(n)
        # Solver
        P = solve_one_phase.pressure(n,P1,P5,k) #pressões tomadas no centróide do elemento

        # Validação
        def P_ans(x): return 1/4*(5-x)
        for x in range(n-1):
            self.assertEqual(np.round(P[x],2),P_ans(x+1),'ValueError')
        plt.figure(0)
        plt.plot(P)
        plt.xlabel('x')
        plt.ylabel('Pressure')
        plt.show()

    def test1D_bifasico(self):
        # Entrada de Dados:
        x0, xf, n, t0, tf, v, CFL = -1., 1., 40, 0., 1., 1., .8

        # Condição Inicial:
        def s(x): return np.sin(np.pi*x)

        #Solver
        vetor = solve_two_phases.saturation(x0,xf,n,s,t0,tf,v,CFL)

    ''' def test_Buckley_Leverett(self):
        x0,xf,n,t0,tf,tv = 0.,1.,40,0.,2.,100
        phi,Area,Qt =.2,1.,0.001# Qt em m3/dia
        CFL = .5
        Swr = 0.16
        Sor = 0.2
        visc_w = .89 #25C
        visc_o = 20*visc_w
        P0 = 1.
        Pl = 0
        v = Qt/Area
        kx = .001 ## talvez tenha uma leve manipulação aqui

        h = 1/n
        t = np.linspace(t0,tf,tv)
        x = np.linspace(x0,xf,n)

        deltaT = CFL*h/v
        A = np.zeros(n)
        B = np.zeros(n)
        C = np.zeros(n)
        T = np.zeros((n,n))
        kro = np.zeros(n)
        krw = np.zeros(n)
        p=0
        Sw = np.zeros(n)
        Sw[0] = 1 #- Sor
        #Sw[1:n] = Swr
        gamma = kx*(deltaT/h**2)

        fig,ax = plt.subplots()
        line, = plt.plot(x, Sw, 'r-')
        ax.set_xlim(x0, xf)
        ax.set_ylim(0, 1)

        while p<tv:
            #Se = (Sw - Swr)/(1 - Swr - Sor)
            krw = Sw**4
            kro = (1-Sw)**2
            #kro[i] = 0.86*math.pow(a,4.8) + 0.04*(1-Se[i])
            #krw[i] = 0.23*math.pow(Se[i],5.5) + 0.18*Se[i]

            lamb_w = krw/visc_w
            lamb_o = kro/visc_o


            for i in range(1,n):
                A[i] = gamma*(lamb_w[i] + lamb_o[i])
                C[i] = gamma*(lamb_w[i-1] + lamb_o[i-1])
                B[i] = A[i] + C[i]
            D = 0

            T[0,0] = 1
            T[n-1,n-1] = 1
            for i in range(1,n-1):
                T[i,i] = B[i]
                T[i,i+1] = -A[i]
                T[i,i-1] = -C[i]
            print('p in',T)
            q = np.zeros(n)
            q[0] = P0
            q[n-1] = Pl

            P = np.matmul(np.linalg.inv(T),q)

            #Sw[0] = S[0] + gamma/phi*(lamb_w[0]*P[0+1] - (lamb_w[0]+lamb_w[n-1])*P[0] + lamb_w[n-1]*P[n-1])
            for i in range(1,n-1):
                Sw[i] = Sw[i] + (gamma/phi)*(lamb_w[i]*P[i+1] - (lamb_w[i]+lamb_w[i-1])*P[i] + lamb_w[i-1]*P[i-1])

            p = p+1

            line.set_ydata(Sw)
            fig.canvas.draw_idle()
            fig.canvas.start_event_loop(0.05)
            plt.show(block=False)
            plt.show(block=False)
        #solve_two_phases.buckley_leverett(x0,xf,n,t0,tf)
'''
