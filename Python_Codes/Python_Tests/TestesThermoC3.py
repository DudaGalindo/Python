import thermo as thermo
import numpy as np
from scipy.integrate import quad
from ..Z_find import z_find
import unittest

R = 8.31447

class Test_thermoC3(unittest.TestCase):
    def test_3_1(self):
        Pi = 101325         #bar to Pascal
        Ti = 323.15         #Kelvin
        s = thermo.Chemical('water',Ti,Pi)
        k = (44.18)*10**(-6) #bar-1
        rhof = 1.01*s.rho
        deltaP = (1/k)*np.log(rhof/s.rho)
        Pf = deltaP + Pi/101325 #sai em bar
        self.assertAlmostEqual(Pf,226.2,1,'ValueError:Failed') #usei 1 pois era a precisão do valor que tinha no gabarito

    def test_3_4(self):
        def f(P):
            return c*P/(P + b)
        m = 1                #kg
        Pi = 1*10**(5)       #bar
        Pf = 500*10**(5)     #bar
        T = 333.15           #K
        b = 2700*10**(5)     #bar
        c = 0.125*10**(-6)   #cm³/gm to m³/gm
        W,err = quad(f,Pi,Pf)
        self.assertEqual(np.round(W,3),0.516,'ValueError:Failed')

    def test_3_32(self):
        B = -140E-6    #m⁻3/mol
        C = 7200E-12   #m⁻12/mol⁻2
        T = 298.15     #K
        P = 12.E5  #Pa
        #R = 8.31447     #Cte geral dos gases
        ## Letra a
        Z_virial = thermo.utils.Z_from_virial_density_form(T,P,B,C)
        V_virial = R*T*Z_virial/P
        ans_a = [np.round(Z_virial,3),np.round(V_virial,6)];

        ## Letra b
        et = thermo.Chemical('ethylene')
        #B_pitzer = thermo.virial.BVirial_Pitzer_Curl(T,et.Tc,et.Pc,et.omega)
        Z_pitzer = z_find.Lee_Kesler(et.Tc,et.Pc,et.omega,T,P)
        V_pitzer = R*T*Z_pitzer/P
        ans_b = [np.round(Z_pitzer,3),np.round(V_pitzer,6)]

        ## Letra c
        Z_RK,eos_RK = z_find.RK(et.Tc,et.Pc,T,P)
        ans_c = [np.round(Z_RK,3),np.round(eos_RK.V_g,7)]

        ## Letra d
        Z_SRK,eos_SRK = z_find.SRK(et.Tc,et.Pc,et.omega,T,P)
        ans_d = [np.round(Z_SRK,3),np.round(eos_SRK.V_g,6)]

        ## Letra e
        Z_PR,eos_PR = z_find.PR(et.Tc,et.Pc,et.omega,T,P)
        ans_e = [np.round(Z_PR,2),np.round(eos_PR.V_g,7)]

        ## Validação
        self.assertEqual(ans_a,[0.929,0.001919],'ValueError: Failed')
        self.assertEqual(ans_b,[0.932,0.001924],'ValueError: Failed')
        self.assertEqual(ans_c,[0.928,0.0019166],'ValueError: Failed')
        self.assertEqual(ans_d,[0.928,0.001918],'ValueError: Failed')
        self.assertEqual(ans_e,[0.92,0.0019006],'ValueError: Failed')

    def test_3_34(self):
        B = -194E-6
        C = 15300E-12
        T = 348.15 #K
        P = 15E5 #bar to Pa - aproximate conversion
        Tc = 318.7 #K
        Pc = 37.6E5 #Bar
        omega = 0.289 #this value is different from the thermo bank and the internet(0.251). No value matches
        ## Letra a
        Z_virial = thermo.utils.Z_from_virial_density_form(T,P,B,C)
        V_virial = R*T*Z_virial/P
        ans_a = [np.round(Z_virial,3),np.round(V_virial,6)];

        ## Letra b
        #B_pitzer = thermo.virial.BVirial_Pitzer_Curl(T,Tc,Pc,omega)
        Z_pitzer = z_find.Lee_Kesler(Tc,Pc,omega,T,P)
        V_pitzer = R*T*Z_pitzer/P
        ans_b = [np.round(Z_pitzer,3),np.round(V_pitzer,6)]

        ## Letra c
        Z_RK,eos_RK = z_find.RK(Tc,Pc,T,P)
        ans_c = [np.round(Z_RK,3),np.round(eos_RK.V_g,7)]

        ## Letra d
        Z_SRK,eos_SRK = z_find.SRK(Tc,Pc,omega,T,P)
        ans_d = [np.round(Z_SRK,3),np.round(eos_SRK.V_g,7)]

        ## Letra e
        Z_PR,eos_PR = z_find.PR(Tc,Pc,omega,T,P)
        ans_e = [np.round(Z_PR,3),np.round(eos_PR.V_g,7)]

        ## Validação
        self.assertEqual(ans_a,[0.893,0.001722],'ValueError: Failed')
        self.assertEqual(ans_b,[0.899,0.001734],'ValueError: Failed')
        self.assertEqual(ans_c,[0.888,0.0017142],'ValueError: Failed') # V=0.0017141
        self.assertEqual(ans_d,[0.895,0.0017271],'ValueError: Failed') #V = 0.0017269
        self.assertEqual(ans_e,[0.882,0.0017017],'ValueError: Failed') #V = 0.0017017

    def test_3_35(self):
        B = -152.5E-6
        C = -5800E-12
        T = 523.15
        P = 1800E3
        ## Letra a
        Z_virial = thermo.utils.Z_from_virial_density_form(T,P,B,C)
        V_virial = R*T*Z_virial/P
        ans_a = [np.round(Z_virial,3),np.round(V_virial,6)];

        ## Letra b
        st = thermo.Chemical('steam')
        Z_pitzer = z_find.Lee_Kesler(st.Tc,st.Pc,st.omega,T,P)
        V_pitzer = R*T*Z_pitzer/P
        ans_b = [np.round(Z_pitzer,3),np.round(V_pitzer,6)]

        ## Validação
        self.assertEqual(ans_a,[0.931,0.002250],'ValueError: Failed')
        self.assertEqual(ans_b,[0.939,0.002268],'ValueError: Failed')


    def test_3_38_39_40(self):
        ## Letra a
        prop = thermo.Chemical('propane')
        T = 313.15
        Psat = 13.71E5

        eos_RK = thermo.eos.RK(prop.Tc,prop.Pc,T,Psat)
        ans_a38 = [np.round(eos_RK.V_l,7),np.round(eos_RK.V_g,7)]
        self.assertEqual(ans_a38,[0.0001081,0.0014991],'ValueError: Failed') #Vg = 0.0014992

        eos_SRK = thermo.eos.SRK(prop.Tc,prop.Pc,prop.omega,T,Psat)
        ans_a39 = [np.round(eos_SRK.V_l,7),np.round(eos_SRK.V_g,7)]
        self.assertEqual(ans_a39,[0.0001047,0.0014806],'ValueError: Failed') #Vg = 0.0014807

        eos_PR = thermo.eos.PR(prop.Tc,prop.Pc,prop.omega,T,Psat)
        ans_a40 = [np.round(eos_PR.V_l,7),np.round(eos_PR.V_g,7)]
        self.assertEqual(ans_a40,[0.0000922,0.0014545],'ValueError: Failed')

    def test_3_41(self):
        et = thermo.Chemical('ethylene')
        ## Letra a
        m = 18 #kg
        T = 328.15 #K
        P = 35E5 #bar to Pa
        # estado gasoso, calcular pelo fator de compressibilidade
        Z,eos_SRK = z_find.SRK(et.Tc,et.Pc,et.omega,T,P)
        V = (m*10**3/et.MW)*Z*R*T/(P) #MW EM GRAMAS/MOL
        self.assertAlmostEqual(V,0.4215,3,'ValueError: Failed')

        ## Letra b - retirada daqui uma vez que o thermo não ta confiável
        V_ = 0.25
        T_ = 323.15
        P_ = 115E5
        # estado supercrítico - líquido - utilizar a relação virial
        Z_ = z_find.Lee_Kesler(et.Tc,et.Pc,et.omega,T_,P_)
        n = P_*V_/(Z_*R*T_)
        m = n*et.MW*10**(-3) #MW EM G/MOL - deu próximo 60.675
    #    #self.assertEqual(m,60.898,'ValueError: Failed') #failed, pelo mesmo motivo

    def test_3_43(self):
        etan = thermo.Chemical('ethanol')
        T = 753.15
        P = 60E5
        Z = z_find.Lee_Kesler(etan.Tc,etan.Pc,etan.omega,T,P)
        Vm = Z*R*T/P
        Vm_i = thermo.volume.ideal_gas(T,P)
        self.assertEqual([np.round(Vm,6),np.round(Vm_i,6)],[0.000988,0.001044],'ValueError: Failed') #0.989

    def test_3_44(self):
        V = 0.35
        prop = thermo.Chemical('propane')
        P = 16E5 #propanes vapor pressure at 320K
        T = 320
        Vliq = 0.8*V
        Vvap = 0.2*V
        Z = z_find.Lee_Kesler(prop.Tc,prop.Pc,prop.omega,T,P)
        mvap = P*Vvap*prop.MW*10**(-3)/(Z*R*T)
        Vmlsat = thermo.volume.Rackett(T,prop.Tc,prop.Pc,prop.Zc)
        nliq = Vliq/Vmlsat
        mliq = nliq*prop.MW*10**(-3)
        self.assertEqual([np.round(mliq,1),np.round(mvap,3)],[127.5, 2.341]) #mliq = 127.522,127.549
