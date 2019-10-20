import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class solve_one_phase1D:
    def permeability(k):
        for i in range(0,len(k)-1):
            if k[i] != k[i+1]:
                if k[i]!=0 and k[i+1]!=0:
                    k[i] = 2*k[i+1]*k[i]/(k[i+1]+k[i]) #ATENÇÃO - CONSIDERANDO QUE A MALHA É UNIFORME
        return k
    def transmissibility(n,k):
        T = np.zeros((n,n))
        T[0,0] = 1
        T[n-1,n-1] = 1
        # Matriz de transmissibilidade --normalizada
        for i in range(1,n-1):
            #t = K[i,i]/h
            T[i,i] = -1*(k[i]+k[i-1])
            T[i,i+1] = 1*k[i]
            T[i,i-1] = 1*k[i-1]
        return T

    def pressure(n,P1,Pn,k):
        q = np.zeros(n)
        q[0] = P1; q[n-1] = Pn  # conhecidos - primeira e quinta linha de [T] determinadas
        k = solve_one_phase1D.permeability(k)
        T = solve_one_phase1D.transmissibility(n,k)
        P = np.matmul(np.linalg.inv(T),q)
        return P

class solve_one_phase2D:
    def permeability(k):
        for i in range(0,len(k)-4):
            if k[i] != k[i+4]:
                if k[i]!=0 and k[i+4]!=0:
                    k[i] = 2*k[i+4]*k[i]/(k[i+4]+k[i]) #ATENÇÃO - CONSIDERANDO QUE A MALHA É UNIFORME
        return k

    def transmissibility(nx,ny,kx,ky):
        T = np.zeros((nx*ny,nx*ny))
        T[0,0] = 1
        T[nx*ny-1,nx*ny-1] = 1
        # Matriz de transmissibilidade --normalizada
        for i in range(1,nx*nx-1):
            T[i,i] = -1*(kx[i+1]+kx[i]) -1*(ky[i+nx]+ky[i])
            T[i,i+1] = 1*kx[i+1]
            T[i,i-1] = 1*kx[i]
            if i < (nx*nx-4): T[i,i+4] = 1*ky[i+nx]
            if i > 4: T[i,i-4] = 1*ky[i]
        return T


    def pressure(nx,ny,P1,Pn,kx,ky):
        q = np.zeros(nx*ny)
        q[0] = P1; q[nx*ny-1] = Pn
        T = solve_one_phase2D.transmissibility(nx,ny,kx,ky)

        P = np.matmul(np.linalg.inv(T),q.T)
        return P

class solve_two_phases:

    def times(t,tf,CFL,h,v):
        tv = -1
        deltaT = CFL*h/v
        while t<tf:
            tv = tv+1
            t = t+deltaT
        return tv

    def saturation(x0,xf,n,s,t0,tf,v,CFL):
        h = (xf - x0)/n # = deltaX
        tv = solve_two_phases.times(t0,tf,CFL,h,v)
        t = np.linspace(t0,tf,tv)
        x = np.linspace(x0,xf,n)
        S = s(x)

        vetor = np.zeros((tv,n)) #inicializing answer vector
        vetor[0,] = S #em t = t0 = 0
        fig,ax = plt.subplots(1)
        line, = plt.plot(x, S, 'r-')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)

        p = 1
        if v>0: #isso é observado pelo gradiente de pressão (se Pi - Pi-1 é positivo)
            while p<tv:
                S[0] = S[0] - CFL*(S[0] - S[n-1]) # The first block is correlated to the last one, so it is calculated separatedly
                for i in range(1,n):
                    S[i] =  S[i] - CFL*(S[i]-S[i-1])
                vetor[p,] = S
                p = p + 1
                line.set_ydata(S)
                fig.canvas.draw_idle()
                fig.canvas.start_event_loop(.05)
                plt.show(block=False)
            #    ani = animation.Func_Animation()
            return vetor

#    def buckley_leverett(x0,xf,n,t0,tf): #TUDO ERRADO
