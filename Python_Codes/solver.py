import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class solve_one_phase1D:
    def conectivity(n):
        conec = np.zeros((n,2));conecy = np.zeros((n,2))
        i = 0
        e = 1
        for el in range(0,n):
            conec[el,0] = e ;conec[el,1] = e+1
            if el == n*(i+1):
                e = e+1
                conec[el,0] = e; conec[el,1] = e+1
                i = i+1
            e = e+1
        conec = conec.astype(int)
        return conec

    def permeability(k,n):
        e = 1
        k_face = np.zeros(len(k)+1)
        for i in range(1,len(k)):
            k_face[e] = k[i]
            if k[i] != k[i-1]:
                k_face[e] = 2*k[i]*k[i-1]/(k[i]+k[i-1])
            e=e+1
            if e==n*(i+1)+i: e=e+2 #ATENÇÃO - CONSIDERANDO QUE A MALHA É UNIFORME
        return k_face

    def transmissibility(n,k,xP): #tentar otimizar usando a matriz esparsa do scipy
        T = np.zeros((n,n))
        conec = solve_one_phase1D.conectivity(n)
        for i in range(0,len(xP)):
            T[xP[i]-1,xP[i]-1] = 1

        # Matriz de transmissibilidade --normalizada
        nl = n-1; init = 1
        if len(xP)==1:
            if xP[0] == 1: init = 1; nl=n
            if xP[0] == n: init = 0

        for i in range(init,nl):
            Face2 = conec[i,1]-1; Face1 = conec[i,0] - 1
            T[i,i] = -1*(k[Face2]+k[Face1])
            if i+1 < n:T[i,i+1] = 1*k[Face2] #quando não tenho CCNewman
            T[i,i-1] = 1*k[Face1]
        return T

    def pressure(n,P1,Pn,xPn,k):
        q = np.zeros(n)
        q[0] = P1; q[n-1] = Pn  # conhecidos - primeira e quinta linha de [T] determinadas
        k = solve_one_phase1D.permeability(k,n)
        T = solve_one_phase1D.transmissibility(n,k,xPn)
        print(T)
        P = np.matmul(np.linalg.inv(T),q)
        return P

    def pressure_Newmann(n,Pn,xPn,qn,xqn,k,L):
        q = np.zeros(n)
        h = L/n
        q[xPn-1] = Pn; q[xqn-1] = qn*h
        k = solve_one_phase1D.permeability(k,n)
        T = solve_one_phase1D.transmissibility(n,k,xPn)
        print('D:',T)
        P = np.matmul(np.linalg.inv(T),q)
        return P

class solve_one_phase2D:
    def conectivity(nx,ny):
        conec = np.zeros((nx*ny,2));conecy = np.zeros((nx*ny,2))
        i = 0; e = 1;x = 0; xold = 0; v = 1
        for el in range(0,nx*ny):
            conec[el,0] = e ;conec[el,1] = e+1
            if el == nx*(i+1):
                e = e+1
                conec[el,0] = e; conec[el,1] = e+1
                i = i+1
            e = e+1

            conecy[el,0] = v ;conecy[el,1] = v+1
            if el == nx*(x):
                x = x+1
                conecy[el,0] = x; conecy[el,1] = x+1

            if x != xold: v = x
            v = nx + v+1
            xold = x

        conec = conec.astype(int)
        conecy = conecy.astype(int)
        return conec,conecy

    def permeability(k,nx,ny):
        k_facex = np.zeros(nx*(nx+1));k_facey = np.zeros(ny*(ny+1))
        e = 1
        for i in range(0,ny):
            for j in range(1,nx):
                k_facex[e] = k[i,j]
                if k[i,j] != k[i,j-1]:
                    k_facex[e] = 2*k[i,j]*k[i,j-1]/(k[i,j]+k[i,j-1])
                e=e+1
                if e==nx*(i+1)+i: e=e+2

        e = 1
        for j in range(0,ny):
            for i in range(1,nx):
                k_facey[e] = k[i,j]
                if k[i,j] != k[i-1,j]:
                    k_facey[e] = 2*k[i,j]*k[i-1,j]/(k[i,j]+k[i-1,j])
                e=e+1
                if e==ny*(j+1)+j: e=e+2

        return k_facex,k_facey

    def transmissibility(nx,ny,kfx,kfy):
        conec,conecy = solve_one_phase2D.conectivity(nx,ny)
        T = np.zeros((nx*ny,nx*ny))
        T[0,0] = 1
        T[nx*ny-1,nx*ny-1] = 1
        # Matriz de transmissibilidade --normalizada
        for el in range(1,nx*ny-1):
            Face2x = conec[el,1]-1; Face1x = conec[el,0]-1
            Face2y = conecy[el,1]-1; Face1y = conecy[el,0]-1
            #Face2y = conec[el]
            T[el,el] = -1*(kfx[Face2x]+kfx[Face1x]) -1*(kfy[Face2y]+kfy[Face1y])
            T[el,el+1] = 1*kfx[Face2x]
            T[el,el-1] = 1*kfx[Face1x]
            if el<nx*ny-nx:T[el,el+4] = 1*kfy[Face2y]
            if el>nx-1:T[el,el-4] = 1*kfy[Face1y]
        return T


    def pressure(nx,ny,P1,Pn,k):
        q = np.zeros(nx*ny)
        q[0] = P1; q[nx*ny-1] = Pn
        kfx,kfy = solve_one_phase2D.permeability(k,nx,ny)
        T = solve_one_phase2D.transmissibility(nx,ny,kfx,kfy)
        P = np.matmul(np.linalg.inv(T),q)
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
