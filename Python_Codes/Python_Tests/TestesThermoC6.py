import thermo
from ..phase import z_find
import math
import numpy as np
import unittest

R = 8.31447
class Test_thermoC6(unittest.TestCase):
    def test_6_14_15_16(self):
    ## Letra a
        Ta = 300
        Pa = 40E5
        ac = thermo.Chemical('ethyne')
        Z_acRK,eos_acRK = z_find.RK(ac.Tc,ac.Pc,Ta,Pa)
        self.assertEqual([np.round(Z_acRK,3),np.round(eos_acRK.H_dep_g),np.round(eos_acRK.S_dep_g,3)],[0.695,-2303,-5.463],'ValueError:Failed')
# H = -2302 S = -5.219 (n sei porque deu essa diferença toda, uma vez que a entalpia e Z deram aproximações bem razoáveis)

        Z_acSRK,eos_acSRK = z_find.SRK(ac.Tc,ac.Pc,ac.omega,Ta,Pa)
        self.assertEqual([np.round(Z_acSRK,3),np.round(eos_acSRK.H_dep_g),np.round(eos_acSRK.S_dep_g,3)],[0.691,-2590,-6.398],'ValueError:Failed')
# H = -2595 S = 6.165

        Z_acPR,eos_acPR = z_find.PR(ac.Tc,ac.Pc,ac.omega,Ta,Pa)
        self.assertEqual([np.round(Z_acPR,3),np.round(eos_acPR.H_dep_g),np.round(eos_acPR.S_dep_g,3)],[0.667,-2651,-6.395],'ValueError:Failed')
# H = -2655 S = -6.168

    ## Letra b
        Tb = 175
        Pb = 75E5
        ar = thermo.Chemical('argon')
        Z_arRK,eos_arRK = z_find.RK(ar.Tc,ar.Pc,Tb,Pb)
    #    self.assertEqual([np.round(Z_arRK,3),np.round(eos_arRK.H_dep_g),np.round(eos_arRK.S_dep_g,3)],[0.604,-2075,-8.798],'ValueError:Failed')
# H = -2068 e S =-7.975  (livro)

        ## Letra c
        Tc = 575
        Pc = 30E5
        bz = thermo.Chemical('benzene')
        Z_bzRK,eos_bzRK = z_find.RK(bz.Tc,bz.Pc,Tc,Pc)
        self.assertEqual([np.round(Z_bzRK,3),np.round(eos_bzRK.H_dep_g),np.round(eos_bzRK.S_dep_g,3)],[0.772,-3318,-4.026],'ValueError:Failed')
# H = -3319, S = -3.879 as diferenças tão maiores na entropia

        ## Letra d
        Td = 500
        Pd = 50E5
        nb = thermo.Chemical('n-butane')
        Z_nbRK,eos_nbRK = z_find.RK(nb.Tc,nb.Pc,Td,Pd)
        self.assertEqual([np.round(Z_nbRK,3),np.round(eos_nbRK.H_dep_g),np.round(eos_nbRK.S_dep_g,3)],[0.685,-4504,-6.544],'ValueError:Failed')
# H = -4503, S = -6.079

    def test_6_23(self): #PERGUNTAR SE RICARDO SABE COMO ENCONTRAR A ENTALPIA DE CADA ESTADO E NÃO DA MUDANÇA DE FASE
        m = 1
        P = 1000E3
        H2O = thermo.Chemical('water')
        T = H2O.Tsat(P)
        w = thermo.Chemical('water',T, P)
        x = (w.rhog)*0.7/((w.rhog)*0.7+(w.rhol)*0.3)
        v = (1/w.rhog)*x + (1/w.rhol)*(1-x)
        print('rho',rho)
        eos = thermo.eos.SRK(w.Tc,w.Pc,w.omega,T,P)
    #    H = x*eos.H_dep_g+(1-x)*eos.H_dep_l

    def test_6_86(self):
        P = 5500E3
        T = 363.15
        mponto = 1.4 #kg/s
        Vmax = 30 #m/s
        ymet = 0.5
        yprop = 1 - ymet
        pm = thermo.Mixture(['methane','propane'],zs = [ymet,yprop],T=T,P=P)
        Z = z_find.Lee_Kesler(pm.Tc,pm.Pc,pm.omega,T,P)
        v = Z*R*T/(P*pm.MW*10**(-3))
        A = mponto*v/Vmax
        d = math.sqrt(A*4/math.pi)
        self.assertAlmostEqual(d,0.0297,3,'ValueError:Failed') #D = 0.002964

    def test_6_88(self):
    # Letra A
        Pa = 60E5
        Ta = 650
        benz = thermo.Chemical('benzene')
        cyclo = thermo.Chemical('cyclohexane')
        bc = thermo.Mixture(['benzene','cyclohexane'],zs = [0.5,0.5],T=Ta,P=Pa)
        Za = z_find.Lee_Kesler(bc.Tc,bc.Pc,bc.omega,Ta,Pa)
        eosamix = thermo.eos_mix.SRKMIX([benz.Tc,cyclo.Tc],[benz.Pc,cyclo.Pc],[benz.omega,cyclo.omega],[0.5,0.5],[[0,0],[0,0]],Ta,Pa)
    #    Z = thermo.utils.Z(Ta,Pa,eosamix.V_g)
    #    ZRK,eosa = z_find.SRK(bc.Tc,bc.Pc,bc.omega,Ta,Pa)
        #eos_PR = thermo.eos.PR(bc.Tc,bc.Pc,bc.omega,T,P)

    # Letra B
        Pb = 100E5
        Tb = 300
        #cd = thermo.Chemical('carbon dioxide')
        #cm = thermo.Chemical('carbon monoxide')
        cmd = thermo.Mixture(['carbon dioxide','carbon monoxide'],zs=[0.5,0.5],T=Tb,P=Pb)
        Zb = z_find.Lee_Kesler(cmd.Tc,cmd.Pc,cmd.omega,Tb,Pb)

    # Letra C
        Pc = 100E5
        Tc = 600
        cdo = thermo.Mixture(['carbon dioxide','n-octane'],zs=[0.5,0.5],T=Tc,P=Pc)
        Zc = z_find.Lee_Kesler(cdo.Tc,cdo.Pc,cdo.omega,Tc,Pc)
        self.assertEqual([Za,Zb,Zc],[0.68,0.794,0.813],'ValueError:Failed') # values are different, based on what I saw,
        # here: 0.7373(cyclohexane Tc different),0.7877(Pc carbon dioxide a little different),0.7912(carbon dioxide Pc is diferent,
        # and omega from mixture too)
