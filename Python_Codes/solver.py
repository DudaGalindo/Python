import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Valido apenas para malhas uniformes
class solve_one_phase1D:
    def conectivity(n):
        conec = np.zeros((n,2));
        for el in range(0,n):
            conec[el,0] = el+1 ;conec[el,1] = el+2
        conec = conec.astype(int)
        return conec

    def permeability(k,n):
        k_face = np.zeros(n+1);

        for i in range(1,n):
            if k[i] != k[i-1]:
                k_face[i] = 2*k[i]*k[i-1]/(k[i]+k[i-1])
            else: k_face[i] = k[i]
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

    def permeability(k,nx,ny):
        k_facex = np.zeros((nx,nx+1));k_facey = np.zeros((ny+1,ny))

        for i in range(0,ny):
            for j in range(1,nx):
                if k[i,j] != k[i,j-1]:
                    k_facex[i,j] = 2*k[i,j]*k[i,j-1]/(k[i,j]+k[i,j-1])
                else: k_facex[i,j] = k[i,j]
                if k[j,i] != k[j-1,i]:
                    k_facey[j,i] = 2*k[j,i]*k[j-1,i]/(k[j,i]+k[j-1,i])
                else: k_facey[j,i] = k[j,i]
        return k_facex,k_facey

    def transmissibility(nx,ny,kfx,kfy):
        conec = solve_one_phase1D.conectivity(nx)
        print(conec)
        T = np.zeros((nx*ny,nx*ny))
        T[0,0] = 1
        T[nx*ny-1,nx*ny-1] = 1
        i=0
        # Matriz de transmissibilidade --normalizada
        for el in range(1,nx*ny-1):
            linha = int(el/4)
            el_linha = el*(np.sign(4-i))**2 - 4*int(el/4)
            Face2x = conec[el_linha,1]-1; Face1x = conec[el_linha,0]-1
            Face2y = conec[linha,1]-1; Face1y = conec[linha,0]-1;
            '''Como el_linha+1 e el_linha coincide com Face2x e Face1x E linha
            coincide com Face2y e Face1y, não precisaria de conec'''
            T[el,el] = -1*(kfx[linha,Face2x]+kfx[linha,Face1x]) -\
                        1*(kfy[Face2y,el_linha]+kfy[Face1y,el_linha])
            T[el,el+1] = 1*kfx[linha,Face2x]
            T[el,el-1] = 1*kfx[linha,Face1x]
            if el<nx*ny-nx:T[el,el+4] = kfy[Face2y,el_linha]
            if el>nx-1:T[el,el-4] = kfy[Face1y,el_linha]


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
