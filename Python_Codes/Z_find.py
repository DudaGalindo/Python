import thermo


class z_find:
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

class flash:
    def Flash(z,K):
        def f(nv): return sum(z*K/(1+nv*(K-1))) - 1
        nv1 = 0.
        nv2 = 1.
        nv = 0.5
        while f(nv)>10**(-10): #4 stages, each stage having 10 search intervals
            dnv = 0.01
            nv1,nv2 = rSolve.root_search(f,nv1,nv2,dnv)
            nv = (nv1 + nv2)/2.
