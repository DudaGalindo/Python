import thermo
import numpy as np
import math
from Python_Codes.rsolve import rSolve

R = 8.31447
class Z_find:
    def Lee_Kesler(Tc,Pc,omega,T,P):
        B = thermo.virial.BVirial_Abbott(T,Tc,Pc,omega) #this virial form its a simple fit to the Lee/Kesner equation
        # The Lee-Kesner equation works for either liquid and gas
        Z = thermo.utils.B_to_Z(B,T,P) #Tc do ciclohaxano ta diferente do que tem no gabarito, muito provavel que seja por isso a diferença do gabarito
        return Z

    def RK(Tc,Pc,T,P):
        eos = thermo.eos.RK(Tc,Pc,T,P)
        Z = thermo.utils.Z(T,P,eos.V_g)
        return Z,eos

    def SRK(Tc,Pc,omega,T,P):
        eos = thermo.eos.SRK(Tc,Pc,omega,T,P)
        Z = thermo.utils.Z(T,P,eos.V_g)
        return Z,eos

    def PR(Tc,Pc,omega,T,P):
        eos = thermo.eos.PR(Tc,Pc,omega,T,P)
        Z = thermo.utils.Z(T,P,eos.V_g)
        return Z,eos

    def Cubic(A,B):
        def f(Z): return Z**3 - (1 - B)*Z**2 + (A - 2*B - 3*B**2)*Z - (A*B - B**2 - B**3)
        def df(Z): return 3*Z**2 - 2*(1 - B)*Z + (A - 2*B - 3*B**2)
        Z1 = 0.
        Z2 = 1. #compressibility factor limits
        Z,i = rSolve.newton_raphson(f,df,Z1,Z2,1e-3,1)
        return Z

class flash:
    def find_v(z,K):
        # iterative proccess to find the number of moles of the vapor phase
        def f(v): return sum(z*K/(1+v*(K-1))) - 1
        v1 = 0.
        v2 = 1.
        v = 0.1
        while np.abs(f(v))>10**(-4): #4 stages, each stage having 10 search intervals
            dv = 0.001
            v1,v2 = rSolve.root_search(f,v1,v2,dv)
            v = (v1 + v2)/2.
        return v #number of moles in the vapor phase

    def y_to_x(y1,A,P1sat,P2sat):
        def f(x1): return x1*(1-1/y1)*np.exp(A*(1-x1)**2)*P1sat + (1-x1)*np.exp(A*(x1)**2)*P2sat
        x11 = 0.
        x12 = 1.
        x1 = 0.5
        while np.abs(f(x1))>10**(-3): #4 stages, each stage having 10 search intervals
            dx1 = 0.001
            x11,x12 = rSolve.root_search(f,x11,x12,dx1)
            x1 = (x11 + x12)/2.
        return x1

class K_value_PR:

    def fb(z,bs):
        b = np.multiply(z,bs)
        return b

    def func_aT(z,aTs,kij):
        aT = 0
        for i in range(3):
            for j in range(3):
                aT = aT + z[i]*z[j]*math.sqrt(aTs[i]*aTs[j])*(1-kij[i][j])
        return aT

    def K_value(Tcs,Pcs,omegas,zs,kij,T,P,K):
        error = [.1,.1,.1]

        while (sum(error) > 0.001):
            n = len(zs)
            Ks = K

            # Setando os vetores como 0
            alphas = np.zeros(n);
            A_linhas = np.zeros((2,n));B_linhas = np.zeros((2,n))
            A = np.zeros(2);B = np.zeros(2)
            b = np.zeros(2);aT = np.zeros(2)
            phi = np.zeros((2,n));K = np.zeros(n)
            Z = np.zeros(2);B_linhas = np.zeros((2,n))
            A_linhas = np.zeros((2,n));z = np.zeros(n)

        ## Cálculo das fases:
            nv,x,y = thermo.activity.flash_inner_loop(zs,Ks)

        ## Calculo dos coeficientes para encontrar Z:
            Trs = np.divide(T,Tcs)
            for i in range(n):
                alphas[i] = (1+(0.37464 + 1.54226*omegas[i]-0.26992*omegas[i]**2)*(1-math.sqrt(Trs[i])))**(2)

            acs = 0.45724*(R**2*np.divide(np.power(Tcs,2),Pcs))
            aTs = np.multiply(acs,alphas)
            bs = 0.07780*R*np.divide(Tcs,Pcs)

            ## Cálculo dos coeficientes para encontrar phi (ambas as fases):

            for p in range(2):
                if p==0: z = y
                else: z = x
                b[p] = sum(K_value_PR.fb(z,bs))
                aT[p] = K_value_PR.func_aT(z,aTs,kij)
                A[p] = aT[p]*P/(R*T)**2
                B[p] = b[p]*P/(R*T)
                Z[p] = Z_find.Cubic(A[p],B[p])
                B_linhas[p,] = bs/b[p]
                for j in range(n):
                    for i in range(n):
                        A_linhas[p,j] = A_linhas[p,j] + 1/aT[p]*(2*math.sqrt(aTs[j]))*(z[i]*aTs[i]**(1/2*(1-kij[i][j])))
                phi[p,] = np.exp(-np.log(Z[p] - B[p]) + (Z[p] - 1)*B_linhas[p,] - A[p]/(2**(1.5)*B[p])*(A_linhas[p,] - B_linhas[p,])*np.log((Z[p] + (math.sqrt(2) + 1)*B[p])/(Z[p] - (math.sqrt(2) - 1)*B[p])))
                K = phi[1,]/phi[0,]

            error = (Ks - K)**2/(Ks*K)
        return K
